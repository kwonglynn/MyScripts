import glob
import ray
from natsort import natsorted
ray.init()

@ray.remote
def count_HILLS(HILLS):
    with open(HILLS, 'r') as fi:
        lines = fi.readlines()
        num = len(lines)
        return (HILLS, num)

HILLS_list = glob.glob("HILLS-*")
HILLS_list = natsorted(HILLS_list)

results_id = []
for HILLS in HILLS_list:
    HILLS_id = ray.put(HILLS)
    result_id = count_HILLS.remote(HILLS_id)
    results_id.append(result_id)

results = ray.get(results_id)
for HILLS, num in results:
    print("{}: {}".format(HILLS, num))

