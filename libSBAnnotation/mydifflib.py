#!/usr/bin/env python
try:
    from Levenshtein import *
    from warnings import warn
    import heapq

    class StringMatcher:
        """A SequenceMatcher-like class built on the top of Levenshtein"""

        def _reset_cache(self):
            self._ratio = self._distance = None
            self._opcodes = self._editops = self._matching_blocks = None

        def __init__(self, isjunk=None, seq1='', seq2=''):
            if isjunk:
                warn("isjunk not NOT implemented, it will be ignored")
            self._str1, self._str2 = seq1, seq2
            self._reset_cache()

        def set_seqs(self, seq1, seq2):
            self._str1, self._str2 = seq1, seq2
            self._reset_cache()

        def set_seq1(self, seq1):
            self._str1 = seq1
            self._reset_cache()

        def set_seq2(self, seq2):
            self._str2 = seq2
            self._reset_cache()

        def get_opcodes(self):
            if not self._opcodes:
                if self._editops:
                    self._opcodes = opcodes(self._editops, self._str1, self._str2)
                else:
                    self._opcodes = opcodes(self._str1, self._str2)
            return self._opcodes

        def get_editops(self):
            if not self._editops:
                if self._opcodes:
                    self._editops = editops(self._opcodes, self._str1, self._str2)
                else:
                    self._editops = editops(self._str1, self._str2)
            return self._editops

        def get_matching_blocks(self):
            if not self._matching_blocks:
                self._matching_blocks = matching_blocks(self.get_opcodes(),
                                                        self._str1, self._str2)
            return self._matching_blocks

        def ratio(self):
            if not self._ratio:
                self._ratio = ratio(self._str1, self._str2)
            return self._ratio

        def quick_ratio(self):
            # This is usually quick enough :o)
            if not self._ratio:
                self._ratio = ratio(self._str1, self._str2)
            return self._ratio

        def real_quick_ratio(self):
            len1, len2 = len(self._str1), len(self._str2)
            return 2.0 * min(len1, len2) / (len1 + len2)

        def distance(self):
            if not self._distance:
                self._distance = distance(self._str1, self._str2)
            return self._distance

    def get_close_matches(word, possibilities, n=3, cutoff=0.6):
        """Copy of get_close_matches from difflib using the StringMatcher instead of the SequenceMatcher"""
        if not n >  0:
            raise ValueError("n must be > 0: %r" % (n,))
        if not 0.0 <= cutoff <= 1.0:
            raise ValueError("cutoff must be in [0.0, 1.0]: %r" % (cutoff,))
        result = []
        s = StringMatcher()
        s.set_seq2(word)
        for x in possibilities:
            s.set_seq1(x)
            if s.real_quick_ratio() >= cutoff and \
               s.quick_ratio() >= cutoff and \
               s.ratio() >= cutoff:
                result.append((s.ratio(), x))

        # Move the best scorers to head of list
        result = heapq.nlargest(n, result)
        # Strip scores for the best n matches
        return [x for score, x in result]
except ImportError:
    # no StringMatcher from python-Levenshtein available
    import difflib
    get_close_matches = difflib.get_close_matches
