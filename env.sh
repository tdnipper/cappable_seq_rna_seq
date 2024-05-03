# env.sh
#!/bin/bash

# Get current directory, user and group IDs for use in scripts

export USER_ID=$(id -u -n)
export GROUP_ID=$(id -g -n)
export BASEDIR=$(dirname "$0")