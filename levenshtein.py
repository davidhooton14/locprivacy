# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 11:46:12 2019

@author: arachnid (https://gist.github.com/Arachnid/491973) and David Hooton
"""

import bisect, numpy as np
from paths import UCToPath

class NFA(object):
  EPSILON = object()
  ANY = object()
  
  def __init__(self, start_state):
    self.transitions = {}
    self.final_states = set()
    self._start_state = start_state
  
  @property
  def start_state(self):
    return frozenset(self._expand(set([self._start_state])))
  
  def add_transition(self, src, input, dest):
    self.transitions.setdefault(src, {}).setdefault(input, set()).add(dest)

  def add_final_state(self, state):
    self.final_states.add(state)
  
  def is_final(self, states):
    return self.final_states.intersection(states)
  
  def _expand(self, states):
    frontier = set(states)
    while frontier:
      state = frontier.pop()
      new_states = self.transitions.get(state, {}).get(NFA.EPSILON, set()).difference(states)
      frontier.update(new_states)
      states.update(new_states)
    return states
  
  def next_state(self, states, input):
    dest_states = set()
    for state in states:
      state_transitions = self.transitions.get(state, {})
      dest_states.update(state_transitions.get(input, []))
      dest_states.update(state_transitions.get(NFA.ANY, []))
    return frozenset(self._expand(dest_states))
  
  def get_inputs(self, states):
    inputs = set()
    for state in states:
      inputs.update(self.transitions.get(state, {}).keys())
    return inputs
  
  def to_dfa(self):
    dfa = DFA(self.start_state)
    frontier = [self.start_state]
    seen = set()
    while frontier:
      current = frontier.pop()
      inputs = self.get_inputs(current)
      for input in inputs:
        if input == NFA.EPSILON: continue
        new_state = self.next_state(current, input)
        if new_state not in seen:
          frontier.append(new_state)
          seen.add(new_state)
          if self.is_final(new_state):
            dfa.add_final_state(new_state)
        if input == NFA.ANY:
          dfa.set_default_transition(current, new_state)
        else:
          dfa.add_transition(current, input, new_state)
    return dfa


class DFA(object):
  def __init__(self, start_state):
    self.start_state = start_state
    self.transitions = {}
    self.defaults = {}
    self.final_states = set()
  
  def add_transition(self, src, input, dest):
    self.transitions.setdefault(src, {})[input] = dest
  
  def set_default_transition(self, src, dest):
    self.defaults[src] = dest
  
  def add_final_state(self, state):
    self.final_states.add(state)

  def is_final(self, state):
    return state in self.final_states
  
  def next_state(self, src, input):
    state_transitions = self.transitions.get(src, {})
    return state_transitions.get(input, self.defaults.get(src, None))

  def next_valid_string(self, input):
    state = self.start_state
    stack = []
    
    # Evaluate the DFA as far as possible
    for i, x in enumerate(input):
      stack.append((input[:i], state, x))
      state = self.next_state(state, x)
      if not state: break
    else:
      stack.append((input[:i+1], state, None))

    if self.is_final(state):
      # Input word is already valid
      return input
    
    # Perform a 'wall following' search for the lexicographically smallest
    # accepting state.
    while stack:
      path, state, x = stack.pop()
      x = self.find_next_edge(state, x)
      if x:
        path += x
        state = self.next_state(state, x)
        if self.is_final(state):
          return path
        stack.append((path, state, None))
    return None

  def find_next_edge(self, s, x):
    if x is None:
      x = u'\0'
    else:
      x = chr(ord(x) + 1)
    state_transitions = self.transitions.get(s, {})
    if x in state_transitions or s in self.defaults:
      return x
    labels = sorted(state_transitions.keys())
    pos = bisect.bisect_left(labels, x)
    if pos < len(labels):
      return labels[pos]
    return None

def levenshtein_automata(term, k): 
  # more complete, slower automata
  nfa = NFA((0, 0))
  for i, c in enumerate(term):
    for e in range(k + 1):
      # Correct character
      nfa.add_transition((i, e), c, (i + 1, e))
      if e < k:
        # Deletion
        nfa.add_transition((i, e), NFA.ANY, (i, e + 1))
        # Insertion
        nfa.add_transition((i, e), NFA.EPSILON, (i + 1, e + 1))
        # Substitution
        nfa.add_transition((i, e), NFA.ANY, (i + 1, e + 1))
  for e in range(k + 1):
    if e < k:
      nfa.add_transition((len(term), e), NFA.ANY, (len(term), e + 1))
    nfa.add_final_state((len(term), e))
  return nfa

def levenshtein_automata2(term, k):
  # removes insertions and deletions but runs faster than above
  nfa = NFA((0, 0))
  for i, c in enumerate(term):
    for e in range(k + 1):
      # Correct character
      nfa.add_transition((i, e), c, (i + 1, e))
      if e < k:
        # Substitution
        nfa.add_transition((i, e), NFA.ANY, (i + 1, e + 1))
  for e in range(k + 1):
    if e < k:
      nfa.add_transition((len(term), e), NFA.ANY, (len(term), e + 1))
    nfa.add_final_state((len(term), e))
  return nfa

def find_all_matches(word, k, lookup_func,*args):
  """Uses lookup_func to find all words within levenshtein distance k of word.
  
  Args:
    word: The word to look up
    k: Maximum edit distance
    lookup_func: A single argument function that returns the first word in the
      database that is greater than or equal to the input argument.
  Yields:
    Every matching word within levenshtein distance k from the database.
  """
  a = args[0] # Q or ll
  # using *args for future reference, nothing else implemented yet
  
  lev = levenshtein_automata2(word, k).to_dfa()
  match = lev.next_valid_string("@")
  while match:
    next = lookup_func(match,a)
    if not next:
      return
    if match == next:
      yield match
      next = next + "@"
    match = lev.next_valid_string(next)
    
    
def lookup_path(word, ll):
    """
    Traverses a tree (implicitly created) of possible words to find the lexicographically 
    smallest valid path greater than or equal to input word (should be a path in UC form)
    """
    
    nrow = len(ll)
    flag = True
    path = [int(x) for x in UCToPath(word).split()]
    
    # find the endpoint index of the longest valid substring
    i = int(path[0] in range(0,nrow))
    while i < len(word) and i > 0:
        if path[i] not in ll[path[i-1]]: break
        i+=1
    
    tmpstr = word[0:i]
    if i == len(word): 
        return(tmpstr)
        
    # increases the character value of the last point by 1
    # if not possible, move back along the string and try again
    while(i >= 0):
        tmpchr = word[i]
        while tmpchr:
            if i > 0:
                # valid characters for word[i] are outgoing edges of word[i-1]
                children = ll[ord(word[i-1])-65]
                children = [chr(x+65) for x in children]
            else: children = [chr(x+65) for x in range(0,nrow)]
            if tmpchr in children and flag:
               tmpstr = tmpstr + tmpchr
               return(tmpstr)
            else:
                tmpchr = nextChar(tmpchr,children)
                flag = True
        i -= 1
        tmpstr = tmpstr[0:-1]
        flag = False
    return False

def nextChar(char,poss):
    """
    given list of unicode characters poss, 'round up' char to smallest p in poss
    such that char <= p, or return False if none exists
    """
    char = ord(char)-65
    poss = [ord(c)-65 for c in poss]
    for i in range(0,len(poss)):
        if char < poss[i]:
            return chr(poss[i]+65)
    return False