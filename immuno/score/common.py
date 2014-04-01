
def memoize(original_fn):
    _cache = {}
    def cached_fn(*args):
        if args in _cache:
            return _cache[args]
        result = original_fn(*args)
        _cache[args] = result
        return result
    return cached_fn