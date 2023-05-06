def flatten(lst):
    result = []
    for item in lst:
        if isinstance(item, list):
            result.extend(flatten(item))
        else:
            result.append(item)
    return result

def list_join(*args, sep=''):
    res = sep.join(flatten(args))
    return res

