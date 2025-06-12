msmd.open('/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a31/group.uid___A001_X1465_X3a32/member.uid___A001_X1465_X3a33/calibrated/calibrated_final.ms')
for spw in msmd.spwsforfield('BrickMaser'):
    frqs = msmd.chanfreqs(spw)
    minf, maxf = frqs.min(), frqs.max()
    print()
    print(f"Spectral Window: {spw}")
    print(f"{minf/1e9}~{maxf/1e9}GHz")
    print()