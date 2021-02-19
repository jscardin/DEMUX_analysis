from pylab import *
def metricinv(v):
    me=['f','p','n','u','m',' ','K','M','G']

    if type(v)==str:
        try:
            return(array(eval(v)))
        except:
            a=[i.isalpha() for i in v]
            if a.count(True):
                out=float(v[:a.index(True)])
                if me.count(v[a.index(True)]):
                    out*=10**(me.index(v[a.index(True)])*3-15)
            else:
                out=float(v)
            return(out)
            
    else:
        return(v)

def metric(i):
    me=['f','p','n','u','m','','K','M','G']
    if i==0:
        a='0'
    elif i==inf:
        a='inf'
    else:
        a=log10(abs(i))
        a=int(a//3)
        a=min(max(a,-5),3)
        b=i*10**(-a*3);
        a='%0.5g'%b+me[a+5]
    return(a)



#print(metricinv('[1 3 4 5 56 6555 4.3 3e6]'))