#ifndef _UTILS_
#define _UTILS_

inline pair<long long,long long> GetStartEnd(int nrank, int rank, long long Min, long long Max)
{
    if(Max < Min) return make_pair(Min, Max);

    long long n = Max - Min + 1;

    //long long start = rank/(nrank+0.) * n;
    //long long end = (rank+1)/(nrank+0.) * n- 1;

    long long start, end;

    if(nrank <= n) {
        start =  (rank*n)   /nrank;
        end   = ((rank+1)*n)/nrank - 1;
    }
    else {
        if(rank < n)
            start = end = rank;
        else {
            start = 0;
            end = -1;
        }

    }
    start += Min;
    end += Min;

    return make_pair(start, end);

}


#endif
