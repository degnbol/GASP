#!/usr/bin/env julia
using ArgParse
using DataFrames, CSV
using Statistics
using Chain: @chain
include("utilities/glob.jl")
include("utilities/AUC.jl")
include("utilities/dataframe.jl")
ROOT = readchomp(`git root`)


parser = ArgParseSettings(description=
"""Train (-o/--train) or evaluate (-i/--eval) a classifier. 
Multiple classifiers are available.""")
@add_arg_table parser begin
    "infiles"
    arg_type = String # default is Any
    nargs = '+'
    help = "Infile with a datapoint on each row for training or prediction,
    the former requiring the presence of a class column (-c/--class).
    Extra files supplied should be feature annotation files,
    which will be (left)joined onto the first (left)file using any shared columns given by -j/--join."
    "--delim", "-d"
    arg_type = Char
    default = '\t'
    "--join", "-j"
    arg_type = String
    nargs = '+'
    help = "Delimiter separating columns in table files."
    "--eval", "-i"
    arg_type = String
    help = "Evaluation: read previously trained model from path."
    "--train", "-o"
    arg_type = String
    help = "Training: train model and write it to path."
    "--class", "-c"
    arg_type = String
    nargs = '+'
    help = "Required for training. Column name containing true class.
    Multiple can be given for getting cor and AUC on testset."
    "--model", "-m"
    arg_type = String
    nargs = '+'
    default = ["randomforest", "n_trees=100"]
    help = "Choose classifier type in case of training.
    First arg is name of classifier, extra args are given to the classifier."
    "--test", "-t"
    help = "For classifiers that support train-test-eval, use this to set test set for early stopping.
    Otherwise it is used for a quick evaluation, i.e. combines train and prediction.
    If a path is given it is assumed to be a table file similar to the first infile.
    If a string is given found as a column name in any of the infiles,
    then the column should contain 0 and 1, where 1 indicates testset,
    or (not implemented yet) more than two unique values which each indicate sets for cross-validation.
    If an integer/fraction is given then this many/fraction of datapoints will be randomly chosen for testing."
    "--substitute", "-s"
    arg_type = String
    nargs = 2
    help = "Encode a column with a substitution matrix.
    First arg is name of column to encode,
    second arg is path to substitution matrix,
    either as absolute path, relative to pwd or to data/NCBI/.
    Should contain numeric columns for substitution and a single string column for naming the new features.
    Default file extension is .tsv
    so simply \"blosum64Amb\" or \"matchAmb\" will be enough."
end

# if run as script
if abspath(PROGRAM_FILE) == @__FILE__
    args = ARGS
else
    # if interactive, do an example
    ROOT = readchomp(`git root`)
    CHEM = glob("$ROOT/results/*-chemicalFeatures") |> only
    pTrain = glob("$ROOT/results/*-features/train.tsv") |> only
    pChems = @chain eglob("$CHEM/{acceptors,20220215}_features.tsv") join(' ')
    pSeqs = glob("$ROOT/results/*-align/muscle_qual.hmm.nterm.tsv") |> only
    pTestset = glob("$ROOT/results/*-experimentalValidation/experimentYield.tsv") |> only
    sArgs = "$pTrain $pChems $pSeqs -j cid enzyme -t $pTestset -c reaction Yield -s seq matchAmb -o hej"
    args = split(sArgs, ' ')
end
args = parse_args(args, parser, as_symbols=true) |> NamedTuple

# check args
@assert (length(args.infiles) < 2) || (args.join !== nothing) "Multiple infiles requires -j/--join"
@assert sum([args.train, args.eval] .!== nothing) == 1 "Provide exactly one of -i/--train or -o/--eval"
@assert (args.train === nothing) || (args.class !== nothing) "-c/--class is required for training"

readdf(source::AbstractString) = CSV.read(source, DataFrame; delim=args.delim)

function dropmissing_verbose!(df, cols=:)
    cols isa Colon || intersect!(cols, names(df))
    colHasMiss = [n for n in names(df[!, cols]) if any(ismissing, df[!, n])]
    before = nrow(df)
    dropmissing!(df, cols)
    after = nrow(df)
    if before > after
        @warn "dropping rows: $before -> $after" colHasMiss
    end
end

dfs = readdf.(args.infiles);
@info "infile(s):" size.(dfs)
df = joinleft(dfs, args.join);
@info "after leftjoins:" size(df)
dropmissing_verbose!(df)

# make sure exactly one of provided class column names is found in the trainset
class = intersect(args.class, names(df)) |> only

if args.test !== nothing
    if isfile(args.test)
        df_test = readdf(args.test);
        @info "test:" size(df_test)
        df_test = joinleft([df_test, dfs[2:end]...], args.join);
        @info "after leftjoins:" size(df_test)
        dropmissing_verbose!(df_test, names(df))
        # Any extra columns in df_test will not be features, and will be 
        # written out unaffected. Extra columns in train should only be the 
        # -c/--class.
        extraInTrain = setdiff(names(df), names(df_test))
        extraInTrainTypes = df[!, extraInTrain] |> eachcol .|> eltype
        # we are only concerned about potential features
        extraInTrain = extraInTrain[extraInTrainTypes .<: Real]
        # class can be missing from test data, we might not know the truth 
        # about something we are predicting on.
        extraInTrain = extraInTrain[extraInTrain .!= class]
        @assert length(extraInTrain) == 0 "Columns missing in -t/--test: $extraInTrain"
    elseif args.test ∈ names(df)
        df_test = df[args.test .== 1]
        df = df[args.test .!= 1]
        select!(df_test, Not(args.test))
        select!(df, Not(args.test))
    else
        try nTest = parse(Int, args.test)
        catch
            try nTest = parse(Float64, args.test) * nrow(df)
            catch
                error("Not understood: -t/--test " + args.test)
            end
        end
        testIdx = shuffle([ones(Int, nTest); zeros(Int, nrow(df) - nTest)])
        df_test = df[testIdx]
        df = df[.!testIdx]
    end
end

"Flatten per sequence."
function substitute(sequence::Vector{Char}, alphabet::Vector{Char}, sub::Matrix)
    dropdims(vcat((sub[:, alphabet .== l] for l in sequence)...); dims=2)
end

# perform substitution
if args.substitute !== nothing
    subCol, subPath = args.substitute
    if !endswith(subPath, ".tsv") subPath *= ".tsv" end
    if !isfile(subPath) subPath = "$ROOT/data/NCBI/" * subPath end
    @assert isfile(subPath) "File not found: $(args.substitute[2])"
    subMat = readdf(subPath)
    # gap is written as both * and -
    if "-" ∉ names(subMat) && "*" ∈ names(subMat) subMat[!, "-"] .= subMat[!, "*"] end
    if "*" ∉ names(subMat) && "-" ∈ names(subMat) subMat[!, "*"] .= subMat[!, "-"] end
    # find single string column
    stringColInd = findall(eltype.(eachcol(subMat)) .<: AbstractString) |> only
    outAlphabet = subMat[!, stringColInd]
    subMat = subMat[!, Not(stringColInd)]
    inAlphabet = names(subMat) .|> only
    subMat = Matrix(subMat)
    
    function substitute(df::DataFrame)
        # make dataframe mapping from unique sequences to their substitutions 
        sequences = unique(df[!, subCol])
        seqlen = length.(sequences) |> unique |> only
        decimals = seqlen |> string |> length
        newSubCols = subCol*"_" .* vcat([outAlphabet .* lpad(i, decimals, '0') for i in 1:seqlen]...)
        subs = hcat(substitute.(collect.(sequences), Ref(inAlphabet), Ref(subMat))...)
        df_sub = DataFrame(subs', newSubCols)
        df_sub[!, subCol] = sequences
        
        # annotate the data with the substitutions. Inner- instead of leftjoin 
        # since it is guranteed that all sequences will be found, since they 
        # are taken from df and this will mean we get Float64 instead of 
        # Float64? (allows missing).
        df = innerjoin(df, df_sub; on=subCol)
        select!(df, Not(subCol))
        @info "after substitutions:" size(df)
        df
    end
    df = substitute(df);
    if args.test !== nothing
        df_test = substitute(df_test);
    end
end

# drop trivial features
trivial = @chain df Matrix all(_ .== _[[1], :]; dims=1) vec
if any(trivial)
    @info "Dropping trivial" names(df)[trivial]
    df = df[!, .!trivial];
end

X, y = df[!, Not([class, args.join...])], df[!, class]
realcols = eltype.(eachcol(X)) .<: Real
if !all(realcols)
    @info "Dropping non-number columns: $(names(X)[.!realcols])"
    X = X[!, realcols];
end;
Xnew = df_test[!, names(X)];


# standardized mats
X = Matrix(X);
Xnew = Matrix(Xnew);
μs = mean(vcat(X, Xnew); dims=1)
σs = std(vcat(X, Xnew); dims=1)
X = @. (X - μs) / σs 
Xnew = @. (Xnew - μs) / σs 

using Flux
using ScikitLearn
using ScikitLearn: fit!
using PyCall
using DecisionTree, Statistics, Random
using ScikitLearn: @sk_import
@sk_import neighbors: KNeighborsClassifier
@sk_import svm: SVC
@sk_import naive_bayes: GaussianNB
@sk_import gaussian_process: GaussianProcessClassifier
@sk_import gaussian_process.kernels: RBF
@sk_import ensemble: GradientBoostingClassifier

classifiers = [
    KNeighborsClassifier(2),
    SVC(kernel="linear", C=0.025),
    SVC(gamma=2, C=1),
    GaussianProcessClassifier(1. * RBF(1.)),
    DecisionTreeClassifier(pruning_purity_threshold=0.8),
    RandomForestClassifier(n_trees=1000),
    AdaBoostStumpClassifier(n_iterations=3),
    GaussianNB(),
    GradientBoostingClassifier(n_estimators=10),
]


for clf in classifiers
    fit!(clf, Matrix(X), y)
    
    pred = try 
        decision_function(clf, Matrix(Xnew))
    catch
        predict_proba(clf, Matrix(Xnew))[:, 2]
    end;
    
    try 
        for class_test in intersect(args.class, names(df_test))
            println("cor = ", cor(pred .== 1, df_test[!, class_test]))
            println("AUC = ", AUC(df_test[!, class_test], pred .== 1))
        end
    catch
    end
    
    for class_test in intersect(args.class, names(df_test))
        println("cor = ", cor(pred, df_test[!, class_test]))
        println("AUC = ", AUC(pred, df_test[!, class_test] .> 50))
    end
end








