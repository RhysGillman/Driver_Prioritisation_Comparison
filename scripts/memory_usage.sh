#!/bin/bash

pid=$1

maxmem=0

while [ -d "/proc/${pid}" ]; do
    mem=`cat /proc/${pid}/status | grep VmRSS | awk '{print $2}'`
    if [[ ${mem} -gt ${maxmem} ]]; then
        maxmem=${mem}
    fi
    sleep 1
done

echo ${maxmem}