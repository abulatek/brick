msmd.open('/blue/adamginsburg/abulatek/brick/symlinks/brick_alma_linesurvey/2019.1.00092.S/science_goal.uid___A001_X1465_X3a41/group.uid___A001_X1465_X3a42/member.uid___A001_X1465_X3a43/calibrated/calibrated_final.ms')
for spw in msmd.spwsforfield('BrickMaser'):
    frqs = msmd.chanfreqs(spw)
    minf, maxf = frqs.min(), frqs.max()
    print()
    print(f"Spectral Window: {spw}")
    print(f"{minf/1e9}~{maxf/1e9}GHz")
    print()