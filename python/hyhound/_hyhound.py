import contextlib, importlib, sys

_variant = "generic"
with contextlib.suppress(ModuleNotFoundError):
    from ._dispatch import get_dispatch_name

    _variant = get_dispatch_name()
_target_name = __name__ + "_" + _variant
_target = importlib.import_module(_target_name, package=__package__)
setattr(_target, "variant", _variant)
sys.modules[__name__] = _target
