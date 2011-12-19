from collections import defaultdict

import math, re, sys

dict_filename = 'n869.txt'

# {{{ digraph counts (anywhere; at word start/end)

print('processing word dictionary {!r}...'.format(dict_filename), file=sys.stderr)
trimmer = re.compile('[^a-zA-Z]')
cnt = defaultdict(lambda: 0)
word_start = defaultdict(lambda: 0)
word_end = defaultdict(lambda: 0)
total_cnt = 0
total_start_end = 0
with open(dict_filename, 'r') as fin:
    for line in fin:
        parts = line.split()
        for part in parts:
            tmp = trimmer.sub('', part).lower()
            l = len(tmp)
            for i in range(l-1):
                cnt[tmp[i:i+2]] += 1
                total_cnt += 1
            if l > 1: # has a digraph at the start and at he end
                word_start[tmp[:2]] += 1
                word_end[tmp[-2:]] += 1
                total_start_end += 1
                
print('done', file=sys.stderr)

# }}}

digraph_prob = {}
digraph_start_prob = {}
digraph_end_prob = {}

# {{{ calculate probabilities

cnt_digraphs = 26*26
for c1 in range(26):
    for c2 in range(26):
        digraph = chr(c1 + 97) + chr(c2 + 97)
        # applying Laplace smoothing
        digraph_prob[digraph] = math.log(cnt[digraph] + 1) - math.log(total_cnt + cnt_digraphs)
        digraph_start_prob[digraph] = math.log(word_start[digraph] + 1) - math.log(total_start_end + cnt_digraphs)
        digraph_end_prob[digraph] = math.log(word_end[digraph] + 1) - math.log(total_start_end + cnt_digraphs)

# }}}

average_digraph_prob = sum(digraph_prob.values()) / cnt_digraphs

stripes_input_filename = 'text_stripes.txt'
orig_text = [] # used for reconstruction at the end
processed_text = [] # used during calculation

# {{{ stripes input processing

print('processing stripes input {!r}...'.format(stripes_input_filename), file=sys.stderr)
with open(stripes_input_filename, 'r') as stripes_fin:
    for line in stripes_fin:
        i = 0
        for element in line.strip('|\n').split('|'):
            if len(processed_text) <= i: # this only happens on the first pass
                processed_text.append([])
                orig_text.append([])
            orig_text[i].append(element)
            tmp = element.lower()
            assert(len(tmp) == 2)
            processed_text[i].append(tmp)
            i += 1
print('done', file=sys.stderr)

# }}}

nstripes = len(orig_text)
if nstripes > 20:
    print('too many stripes ({}) for dynamic programming'.format(nstripes), file=sys.stderr)
    sys.exit(-1)
stripe_len = len(orig_text[0])

# dynamic programming solution {{{

def join_digraphs(a, b):
    if a[1] == ' ':
        if a[0] == ' ': # two spaces are unlikely before anything
            if b[0] == ' ':
                return average_digraph_prob - 5 
            else: # especially before a nonspace
                return average_digraph_prob - 10 
        elif b in digraph_start_prob:
            return digraph_start_prob[b]
        else:
            return average_digraph_prob 
    elif b[0] == ' ':
        if a in digraph_end_prob:
            return digraph_end_prob[a]
        else:
            return average_digraph_prob
    else:
        tmp = a[1] + b[0]
        if tmp in digraph_prob:
            return digraph_prob[tmp]
        else:
            return average_digraph_prob

join_memo = [[None]*nstripes for _ in range(nstripes)]
def join_stripes(aind, bind):
    global nstripes, stripe_len, join_memo, processed_text
    if join_memo[aind][bind] is None:
        a = processed_text[aind]
        b = processed_text[bind]
        p = 0.0
        for i in range(stripe_len):
            p += join_digraphs(a[i], b[i])
        join_memo[aind][bind] = p
    return join_memo[aind][bind]

memo = [[None]*(1<<nstripes) for _ in range(nstripes)]
next_choice = [[None]*(1<<nstripes) for _ in range(nstripes)]
def get_prob(at, unused):
    global nstripes, stripe_len, memo
    if unused == 0:
        return 1.0
    if memo[at][unused] is None:
        best = -1e128
        chosen_next = -1
        for next in range(nstripes):
            if unused & (1<<next):
                val = join_stripes(at, next) + get_prob(next, unused ^ (1<<next))
                if val > best:
                    best = val
                    chosen_next = next
        assert(chosen_next != -1)
        memo[at][unused] = best
        next_choice[at][unused] = chosen_next
    return memo[at][unused]

# }}}

print('calculating most probable order', file=sys.stderr)
best = -1e128
at = -1
for i in range(nstripes):
    print('{}/{}'.format(i+1, nstripes), file=sys.stderr)
    val = get_prob(i, ((1<<nstripes)-1) ^ (1<<i))
    if val > best:
        best = val
        at = i
print('done', file=sys.stderr)

order = [at]
unused = ((1<<nstripes)-1) ^ (1<<at)
while next_choice[at][unused] is not None:
    next = next_choice[at][unused]
    order.append(next)
    assert(unused & (1<<next))
    at, unused = next, (unused ^ (1<<next))

for row in range(stripe_len):
    print(''.join(orig_text[order[i]][row] for i in range(nstripes)))

# Shannon published his paper in 1948.
