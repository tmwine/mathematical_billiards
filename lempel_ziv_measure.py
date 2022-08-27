'''

given a file of symbolic sequences, assigns a value to each according to
their LZ77 complexity; LZ77 here is a minor variant of the original
Lempel-Ziv paper; complexity is estimated in bytes to compress, using a fixed
number of bits for the lookback index, and a fixed number of bits for the
subsequence length

this requires specification of 2 files in header area:
        in_file--a file of symbolic sequences, one on each line
        out_file--a file to save the LZ77 measure of each
bounce sequence in, one entry per line

note, out_file actually holds the quantity c-LZ77 on each line, where c is
the maximum complexity over all sequences--this is so the most
compressed sequence will have the highest value

also, for optimal compression estimates, especially for periodic sequence
with a long period, set SYM_STR_LEN in header area to length of symbolic
sequences to be processed

'''

import math


#@@@ set sequence length and I/O filenames:
SYM_STR_LEN = 500   # for estimation of efficient bytes to compress, assuming
# full past-search / lookahead windows are wanted
in_file = ""
out_file = ""
#@@@


if in_file=="" or out_file=="":
    print("enter names of both input and output files in script header")
    raise ValueError


SB_LEN = SYM_STR_LEN # search buffer length; how far back from current position
# we're allowed to look to find a substring match
LA_LEN = SYM_STR_LEN # lookahead buffer length; max distince ahead from the
# current position we can look for substrings already covered

# alterantively to allowing computation from SYM_STR_LEN, can just
# ~manually set the buffer lengths
#SB_LEN = 31
#LA_LEN = 31

# how many bits, max, required to record past offset index, and subsequence
# length index:
bits_sb = math.ceil(math.log(SB_LEN)/math.log(2))
bits_la = math.ceil(math.log(LA_LEN)/math.log(2))

def stor_func(tup):
    if tup[0]==0:
        return bits_sb+bits_la+8
    else:
        return bits_sb+bits_la

def seq_cpr(in_st):
    # LZ77 compression function; returns a list of tuples:
        # (lookback steps,subsequence length,first character)
    out_tps = []
    ii = 0
    while ii < len(in_st):
        sbl_now = min(ii,SB_LEN) # how far back are we allowed to look from cursor
        # ii position
        lal_now = min(len(in_st)-ii,LA_LEN) # how far forward are we allowed
        # to look from cursor ii position
        best_match = None
        cur_chr = in_st[ii]
        pst_fnd = True
        for jj in [x+1 for x in range(lal_now)]:
            search_buffer = in_st[ii-sbl_now:ii+jj]
            lookahead_buffer = in_st[ii:ii+jj]
            mtch_indx = search_buffer.find(lookahead_buffer)
            if mtch_indx == -1 or mtch_indx > sbl_now:
                break # substring is not found before or at cursor position (ii)
            elif mtch_indx < sbl_now:
                best_match = (sbl_now-mtch_indx,jj) # substring is found, starting
                # before cursor position (ii)
            elif jj==1:
                pst_fnd = False # if character at cursor position does not have
                # previous matches; preps for special case of new, repeated letter
            else:   # substring match found exactly at cursor position
                if pst_fnd == False and lookahead_buffer == cur_chr*jj:
                    best_match = (0,jj) # handles case of new, repeated letter
                    continue
                else:
                    break
        else:
            jj += 1     # kludge, to advance ii by the correct amount
        if best_match is None:
            out_rng = (0,1,cur_chr)
            ii += 1
        else:
            out_rng = (*best_match,cur_chr)
            ii += jj-1
        out_tps.append(out_rng)
    return out_tps

val_vec = []
with open(in_file) as fp:
    for line in fp:
        out_tps = seq_cpr(line.strip())
        out_bytes = math.ceil(sum([stor_func(tup) for tup in out_tps])/8.0)
        val_vec.append(out_bytes)

mx_vl = max(val_vec)
with open(out_file,'w') as fp:
    for val in val_vec:
        tmp = fp.write(str(mx_vl-val)+'\n')
