#!/usr/bin/python

import sys

from math import log

import qcmdpc_stats
import ci


ALPHA = 0.01
ALGO_PARAM = {"BP": ['bp_scale', 'bp_saturate'],
              "GRAY_B": ['threshold_c0', 'threshold_c1'],
              "GRAY_BGF": ['threshold_c0', 'threshold_c1'],
              "GRAY_BGB": ['threshold_c0', 'threshold_c1'],
              "GRAY_BG": ['threshold_c0', 'threshold_c1'],
              "BACKFLIP2": ['threshold_a0', 'threshold_a1', 'threshold_a2',
                            'threshold_a3', 'threshold_a4', 'ttl_saturate'],
              "BACKFLIP": ['ttl_c0', 'ttl_c1', 'ttl_saturate'],
              "SORT": ['gray_size']}


if len(sys.argv) <= 3:
    print("{} filename1 filename2 block_length3 [alpha]".format(sys.argv[0]))
    sys.exit(1)

filename1 = sys.argv[1]
filename2 = sys.argv[2]
block_length3 = int(sys.argv[3])
alpha = float(sys.argv[4]) if len(sys.argv) > 4 else ALPHA

data1 = qcmdpc_stats.get_data(filename1)
data2 = qcmdpc_stats.get_data(filename2)

if data2['block_length'] < data1['block_length']:
    data1, data2 = data2, data1


for i, data in enumerate([data1, data2]):
    print("Point {}".format(i + 1))
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

    print()


print("Extrapolated values:")
print("{:13}: {}".format('block_length', block_length3))

density = None
distance = None
if all([data1[p] == data2[p] for p in ['index',
                                       'block_weight',
                                       'error_weight',
                                       'weak',
                                       'weak_p',
                                       'error_floor',
                                       'error_floor_p']]):
    n0, d, t, weak, error_floor, m, l = (data1[p]
                                         for p in ['index',
                                                   'block_weight',
                                                   'error_weight',
                                                   'weak',
                                                   'error_floor',
                                                   'weak_p',
                                                   'error_floor_p'])
    density, distance = qcmdpc_stats.get_density(n0,
                                                 block_length3,
                                                 d,
                                                 t,
                                                 weak,
                                                 error_floor,
                                                 m,
                                                 l)

if density != 0:
    print("{:13}: {}".format('density', density))
if distance:
    print("{:13}: {}".format('distance', distance))

cumul1, total1 = qcmdpc_stats.cumul_iteration(data1['iteration'])
cumul2, total2 = qcmdpc_stats.cumul_iteration(data2['iteration'])

min_it = max(min(cumul1.keys()),
             min(cumul2.keys()))
max_it = min(max(cumul1.keys()),
             max(cumul2.keys()))

block_length1 = data1['block_length']
block_length2 = data2['block_length']

dfr = {}
print("{:13}:".format('dfr (with CI)'))
for it in range(min_it, max_it + 1):
    _, f1 = max([(i, f) for (i, f) in cumul1.items() if i <= it])
    _, f2 = max([(i, f) for (i, f) in cumul2.items() if i <= it])
    dfr[it] = tuple(map(lambda x: x / log(2), ci.conf_int(block_length3,
                                                          block_length1, f1, total1,
                                                          block_length2, f2, total2,
                                                          alpha)))
    print("{:16}: {:.3f} {:.3f} {:.3f}".format(it, *dfr[it]))

if density:
    print("{:13}:".format('dfr+density (with CI)'))
    for it in range(min_it, max_it + 1):
        ddfr = tuple(map(lambda x: x + density, dfr[it]))
        print("{:16}: {:.3f} {:.3f} {:.3f}".format(it, *ddfr))
