#!/usr/bin/python

import sys

import qcmdpc_stats


ALPHA = 0.01
ALGO_PARAM = {"BP": ['bp_scale', 'bp_saturate'],
              "GRAY_B": ['threshold_c0', 'threshold_c1'],
              "GRAY_BGF": ['threshold_c0', 'threshold_c1'],
              "GRAY_BGB": ['threshold_c0', 'threshold_c1'],
              "GRAY_BG": ['threshold_c0', 'threshold_c1'],
              "BACKFLIP2": ['threshold_a0', 'threshold_a1', 'threshold_a2', 'threshold_a3', 'threshold_a4', 'ttl_saturate'],
              "BACKFLIP": ['ttl_c0', 'ttl_c1', 'ttl_saturate'],
              "SORT": ['gray_size']}


if len(sys.argv) <= 1:
    print("{} filename [alpha]".format(sys.argv[0]))
    sys.exit(1)

filename = sys.argv[1]
alpha = float(sys.argv[2]) if len(sys.argv) > 2 else ALPHA
data = qcmdpc_stats.get_data(filename)
data['dfr'] = qcmdpc_stats.get_dfr(data['iteration'], alpha)

for p in ['index', 'block_length', 'block_weight', 'error_weight', 'algo']:
    print("{:13}: {}".format(p, data[p]))

if data['algo'] in ALGO_PARAM:
    for p in ALGO_PARAM[data['algo']]:
        print("{:13}: {}".format(p, data[p]))

if data['weak'] != 0:
    for p in ['weak', 'weak_p']:
        print("{:13}: {}".format(p, data[p]))

if data['error_floor'] != 0:
    for p in ['error_floor', 'error_floor_p']:
        print("{:13}: {}".format(p, data[p]))

if data['density'] != 0:
    print("{:13}: {}".format('density', data['density']))

if data['distance']:
    print("{:13}: {}".format('distance', data['distance']))

print("{:13}:".format('dfr (with CI)'))
for it, dfr in sorted((data['dfr'].items())):
    print("{:16}: {:.3f} {:.3f} {:.3f}".format(it, *dfr))

if data['density'] != 0:
    print("{:13}:".format('dfr+density (with CI)'))
    for it, dfr in sorted((data['dfr'].items())):
        print("{:16}: {:.3f} {:.3f} {:.3f}".format(
            it, *list(map(lambda x: x + data['density'], dfr))))
