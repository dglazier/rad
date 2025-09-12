#!/bin/bash

source "/home/garyp/eic/setup_rad.sh"
cd "/home/garyp/eic/rad/examples/"

root -l -b -q "DDVCS_GenHeli.C(\"/w/work5/home/garyp/18x275_ddvcs_1M_events_plus.root\")"
root -l -b -q "DDVCS_GenHeli.C(\"/w/work5/home/garyp/18x275_ddvcs_1M_events_minus.root\")"

mv "/w/work5/home/garyp*_flat.root" "/w/work5/home/garyp/eic/Farm/data/EpIC_ep_DDVCS_18x275/"
