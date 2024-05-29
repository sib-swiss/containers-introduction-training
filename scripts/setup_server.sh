#!/usr/bin/env bash

./generate_credentials -l ../participant_list.txt -o credentials

# taking ubuntu 22.04 image
# ubuntu 24.04 image does not allow for the right permissions apptainer
./multi_instance \
-o credentials \
-t c5a.4xlarge \
-a ami-026c3177c9bd54288 \
-s sg-0b638dae2ff2643d2 \
-b 500 \
-k key_ubuntu_gvg \
-p $SIBKEY 
