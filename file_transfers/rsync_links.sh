#!/usr/bin/env bash

# Script name: rsync_link.sh
# Description:
#
#   Transfers a directory with many symlinks to a new path while
#   retaining the relative relationships between link and source.
#   Source files are also copied and links point to them.
#
#   Note: this will replace all fullpath links with relative ones in
#   your source directory before rsyncing to the new location.
#
# Author: Ryan Kuster
# Date: 2025-02-11
# Usage:
#
#   rsynch_link.sh <source directory> <new directory name>
#
#   If a source file links to a origin file outside of source directory,
#   the file will be copied and the link will not be preserved in the
#   new directory
#
#   The new directory cannot be existing before running.

directory_a=$1
directory_b=$2


# Check if both arguments are provided and not the same.

if [[ -z $directory_a || -z directory_b ]]; then
    echo "Error: Both arguments must be provided."
    exit 1
elif [[ $directory_a == $directory_b ]]; then
    echo "Error: Arguments must not be the same."
    exit 1
elif [[ ! -d $directory_a ]]; then
    echo "Error: $directory_a is not a valid directory."
    exit 1
elif [[ -d $directory_b ]]; then
    echo "Error: $directory_b exists"
    exit 1
elif [[ "$(realpath $directory_a)" == "$(realpath $directory_b)" ]]; then
    echo "Error: Directory paths must not be the same."
    exit 1
fi

echo "Valid inputs: '$1' and '$2'"


# Recursively search a directory (directory_a) for symlinks.
# If a symlink is an absolute path, reassign the symlink to be a
# relative path.

find $directory_a -type l | while read link; do
    target=$(readlink "$link")

    # Check if the symlink is an absolute path.
    if [[ "$target" = /* ]]; then
        rel_target=$(realpath --relative-to="$(dirname "$link")" "$target")
        ln -snf "$rel_target" "$link"
    fi
done

rsync -a --copy-unsafe-links $directory_a $directory_b
