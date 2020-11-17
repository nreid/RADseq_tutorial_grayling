#!/bin/bash

jid1=$(sbatch --parsable 00_symlinkdata.sh )
jid2=$(sbatch --parsable --dependency=afterok:$jid1 01_process_radtags.sh )
