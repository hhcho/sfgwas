#!/bin/sh
ps -a | grep sfgwas | awk '{print $1}' | xargs -I job kill job
