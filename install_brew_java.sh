#!/usr/bin/env zsh
brew tap AdoptOpenJDK/openjdk
brew install --cask adoptopenjdk8
# as per `brew info openjdk` we symlink to make runtime accessible
sudo ln -sfn /usr/local/opt/openjdk/libexec/openjdk.jdk /Library/Java/JavaVirtualMachines/openjdk.jdk
