#!/usr/bin/env bash

timedatectl set-timezone America/New_York
timedatectl status

chown -R ubuntu /home/ubuntu

echo "ubuntu:2019neuro!!" | chpasswd  # update default password!
history -c

# log file: /var/log/cloud-init-output.log
