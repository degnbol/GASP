#!/usr/bin/env python3
import subprocess

def run(cmd):
    """
    Run shell command.
    :param cmd: string command to run
    :return: string stdout of command (without final newline).
    """
    return subprocess.run(cmd.split(), capture_output=True, text=True).stdout.strip()


def git_root(relpath=""):
    """
    Get absolute path given a path relative to git root.
    :param relpath: string path relative to git root.
    :return: string absolute path.
    """
    return os.path.join(run("git root"), relpath)

