import sys

import pandas as pd



# check size of a object
def get_size_mb(obj):
    """Returns the shallow size of an object in megabytes."""
    size_bytes = sys.getsizeof(obj)
    return size_bytes / (1024 ** 2)


# check ditcts if equal:
def dicts_equal(d1, d2):
    if d1.keys() != d2.keys():
        return False
    for key in d1:
        v1, v2 = d1[key], d2[key]
        if isinstance(v1, pd.DataFrame) and isinstance(v2, pd.DataFrame):
            if not v1.equals(v2):
                return False
        else:
            if v1 != v2:
                return False
    return True
