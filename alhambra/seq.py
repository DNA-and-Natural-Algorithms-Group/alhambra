# Utility functions for sequences
import string

_L_TO_N = { 'a': frozenset((0,)),
            'b': frozenset((1, 2, 3)),
            'c': frozenset((1,)),
            'd': frozenset((0, 2, 3)),
            'g': frozenset((2,)),
            'h': frozenset((0, 1, 3)),
            'k': frozenset((2, 3)),
            'm': frozenset((0, 1)),
            'n': frozenset((0, 1, 2, 3)),
            's': frozenset((1, 2)),
            't': frozenset((3,)),
            'v': frozenset((0, 1, 2)),
            'w': frozenset((0, 3)) }

_N_TO_L = { v: i for i,v in _L_TO_N.items() }

def is_null(seq):
    """Return True if a sequence consists only of Ns, or is empty. 
    Return False otherwise."""
    return set(seq.lower()).issubset(set('n').union(set(string.whitespace)))

def is_definite(seq):
    """Return True if a sequence consists only of defined bases.  Return
    False otherwise.  If blank, return False.
    """
    if set(seq.lower()).issubset(set(string.whitespace)):
        return False
    return set(seq.lower()).issubset({'a','g','c','t'}.union(set(string.whitespace))) 

def merge(seq1, seq2):
    """Merge two sequences together, returning a single sequence that
    represents the constraint of both sequences.  If the sequences
    can't be merged, raise a MergeConflictError.

    FIXME: this needs to intelligently handle case and whitespace.
    """
    
    if len(seq1) != len(seq2):
        raise MergeConflictError(seq1,seq2,'length',len(seq1),len(seq2))
    
    out = []
    for i,(n1,n2) in enumerate(zip(seq1.lower(),seq2.lower())):
        try:
            out.append(_N_TO_L[frozenset.intersection(_L_TO_N[n1], _L_TO_N[n2])])
        except KeyError as e:
            if e.args[0] == frozenset():
                raise MergeConflictError(seq1,seq2,i,n1,n2) from None
            else:
                raise e
    return ''.join(out)


class MergeConflictError(ValueError):
    """
    Merge of items failed because of conflicting information.
    Arguments are (item1, item2, location or property, value1, value2)
    """
    
class MergeConflictsError(ValueError):
    """
    Merge of multiple items failed because individual merges
    raised MergeConflictErrors.
    Arguments are ([list of MergeConflictErrors])
    """
