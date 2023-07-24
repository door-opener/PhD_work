import os, sys
import numpy as np
import subprocess as sp

def monthToNum(shortMonth):
    return {
            'January': 1,
            'February': 2,
            'March': 3,
            'April': 4,
            'May': 5,
            'June': 6,
            'July': 7,
            'August': 8,
            'September': 9, 
            'October': 10,
            'November': 11,
            'December': 12
    }[shortMonth]

def readData( in_file ):
    names, code, ID = [], [], []
    f=open(in_file, "r")
    next(f)
    for line in f:
        split = line.split(",")
        names.append( [ split[3].strip(), split[4].strip() ] )
        code.append( split[1].strip() )
        ID.append( split[-1].strip() )
    f.close()
    return names, code, ID

def getURL( url ):
    err = sp.check_output( ('lynx -width=200 -dump "%s"'%(url)), stderr=sp.STDOUT, shell=True).decode('cp850')
    return err

def getGrantURLs( page ):
    page = page.split()
    grant_IDs, grant_URLs = [], []
    for i in range(len(page)):
        if "/" in page[i] and page[i][0] == "[":
           grant_IDs.append( page[i].split("]")[-1])
    for i in range(len(page)):
        if "https://gow.epsrc.ukri.org/NGBOViewGrant" in page[i]:
           for j in range(len(grant_IDs)):
               if grant_IDs[j] in page[i]:
                  grant_URLs.append( page[i] )
    return grant_IDs, grant_URLs

def getGrantTitle( grant_URL ):
    dump = getURL( grant_URL ).split()[0:50]
    rst, cnt, out = False, 0, []
    for i in range(len( dump )):
        if "Title:" in dump[i]:
           rst = True
        if rst:
           if "Principal" in dump[i]:
              break
           out.append( dump[i] )
    return out[1:-1]

def getStEndValue( grant_URL ):
    dump = getURL( grant_URL )
    dump, out = dump.split()[0:190], []
    rst, cnt = False, 0
    for i in range(len(dump) ):
        if "Starts:" in dump[i]:
           rst = True
        if rst:
           cnt +=1
           out.append( dump[i] )
        if cnt == 11:
           break
    return out

def main():

    try:
      in_file = sys.argv[1]
    except:
      print("\nEnter Input FILE!\n")
      sys.exit()

    names, code, ID = readData( in_file )

    for i in range(len(names)):
        if ID[i]:
           print("---------------------------------------------------------------------------------")
           print("\n%s. %s %s" %(names[i][0][0], names[i][1], ID[i]) )
           print("---------------------------------------------------------------------------------")

           dump = getURL( ID[i] )
           grant_IDs, grant_URLs = getGrantURLs( dump )
           sum_value = 0

           for j in range(len( grant_URLs )):
               check, title = getStEndValue( grant_URLs[j] ), getGrantTitle( grant_URLs[j] )
               print(check, title)
               stdate, end_date, value = [ check[1], check[2], check[3] ], [ check[5], check[6], check[7] ], check[-1].replace(",","")

               if monthToNum(stdate[1]) >= 11 and int(stdate[2]) >= 2018:
                  print("\n%s: %s" %(grant_IDs[j], " ".join(title)) )
                  print("%s - %s %s\n" %(" ".join(stdate), " ".join(end_date), check[-1]) )
                  sum_value += int(value)

           print("---------------------------------------------------------------------------------")
           print("Total Value: %d.00"%(sum_value))
           print("---------------------------------------------------------------------------------\n")

if __name__ == "__main__":
   main()
