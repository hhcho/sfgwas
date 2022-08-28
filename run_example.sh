#!/bin/sh
NUM_MAIN_PARTY=2
LOG_PREFIX=stdout
for (( i = 0; i <= $NUM_MAIN_PARTY; i++ )) 
do
  echo "Running PID=$i"
  CMD="PID=$i go run sfgwas.go | tee /dev/tty > ${LOG_PREFIX}_party${i}.txt"
  if [ $i = $NUM_MAIN_PARTY ]; then
    eval $CMD
  else
    eval $CMD &
  fi
done