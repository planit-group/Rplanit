#!/bin/bash

#source @gate.setenv
Gate mac/main.mac

# sanitizza tutto...
for i in output/*.hdr
do sanitize_hdr "$i" "$i"; done
