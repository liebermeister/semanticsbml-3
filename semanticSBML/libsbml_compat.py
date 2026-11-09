"""
Compatibility of the 2010 code with the current libSBML python bindings.

The original code uses libSBML objects (species, reactions ...) as keys of
dictionaries and as members of sets. Under Python 2 this worked; under
Python 3 the libSBML classes define equality (two python objects are equal
if they wrap the same underlying libSBML object) but no hash value, so
they cannot be used as keys any more.

Here the missing hash function is added. It is derived from the address of
the wrapped libSBML object and therefore consistent with libSBML's equality.
"""

import inspect
import libsbml


def _hash(self):
    try:
        return hash(int(self.this))
    except Exception:
        return id(self)


def install():
    for name, cls in inspect.getmembers(libsbml, inspect.isclass):
        if '__eq__' in vars(cls) and vars(cls).get('__hash__', None) is None:
            try:
                cls.__hash__ = _hash
            except TypeError:
                pass


install()


def set_id(element, new_id):
    """Set the identifier of a libSBML element.

    For rules and initial assignments the 'identifier' is the variable they
    define. Old libSBML versions changed that variable in setId(); current
    versions refuse (and leave the element unchanged), so it is done here."""
    if isinstance(element, libsbml.Rule):
        if isinstance(element, libsbml.AlgebraicRule):
            return 0
        return element.setVariable(new_id)
    if isinstance(element, libsbml.InitialAssignment):
        return element.setSymbol(new_id)
    return element.setId(new_id)
