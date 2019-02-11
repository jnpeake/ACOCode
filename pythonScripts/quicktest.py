def get_stats( arr ):
    num = len(arr)
    mean = 0.0
    min = 1e20
    max = -1e20
    for v in arr:
        mean += v
        if v > max:
            max = v
        if v < min:
            min = v
            
    mean /= float(num)
    lq = sorted(arr)[num//4]
    median = sorted(arr)[num//2]
    uq = sorted(arr)[3*num//4]
    perAnt = mean/1024
    print ("Num: ",num, " Min: ",min, " LQ: " ,lq," Median: ",median, " UQ: ",uq," Max: ",max," Mean: ",mean)
    return (min, lq, median, uq, max, mean)


