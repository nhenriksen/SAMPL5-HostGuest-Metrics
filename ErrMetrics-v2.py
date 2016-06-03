import sys,re
import numpy as np
from scipy import stats

OAHOnly = False
OAMeOnly = False
CalcPairDiffs = False
CorrectOA = False
CorrectCB = False
WithUncert = True
WithRep = True

### Arguments: <Flags> <ExpType> <ExperimentFile> <CalculationFile1> [<CalculationFile2> ...]
### Note, files should have the same number of data points
calcfiles=[]
for arg in sys.argv:
  if arg == 'OAHOnly':
    OAHOnly = True
  if arg == 'OAMeOnly':
    OAMeOnly = True
  if arg == 'CalcPairDiffs':
    CalcPairDiffs = True
  if arg == 'CorrectOA':
    CorrectOA = True
  if arg == 'CorrectCB':
    CorrectCB = True
  if re.search(r'^(ka|dg|dh)$', arg):
    exptype = arg
  if re.search(r'\.txt$', arg):
    if re.search(r'Exp', arg):
      expfile = arg
    else:
      calcfiles.append(arg)

### Clean up method names
CalcNames=[]
for name in calcfiles:
  cols1 = name.split('.')
  cols2 = cols1[0].split('/')
  CalcNames.append(cols2[-1])
  CalcNames[-1] = re.sub('TIxBAR','TI/BAR', CalcNames[-1])

### Settings
R = 0.0019872036
T = 298.0
BootCyc = 100000
MNames = ('Slope', 'Interc', 'R', 'R^2', 'RMSE', 'MSE', 'MUE', 'TAU')
Nm = len(MNames)

### Error Metric Calculations
def geterrormetrics(x, y):
  MTmp = np.zeros([Nm], np.float64)
  # abs(1-Slope), Intercept, R
  Slp, MTmp[1], MTmp[2], pval, stderr = stats.linregress(x,y)
  MTmp[0] = Slp 
  #MTmp[0] = np.abs(1-Slp)
  # R^2
  MTmp[3] = MTmp[2]**2
  # RMSE
  MTmp[4] = np.sqrt(np.mean(((y - x)**2)))
  # MSE
  MTmp[5] = np.mean((y - x))
  # MUE
  MTmp[6] = np.mean(np.absolute(y - x))
  # Tau
  MTmp[7], prob = stats.kendalltau(x,y)
  if np.isnan(MTmp[7]):
    MTmp[7] = 0.0

  return (MTmp)


### BootStrapping Definition
def bootstrap(x, xsem, y, ysem):
  MBoot = np.zeros([Nm,BootCyc], np.float64)
  MVals = np.zeros([Nm], np.float64)
  MSEMs = np.zeros([Nm], np.float64)
  xtmp = np.zeros([len(x)], np.float64)
  ytmp = np.zeros([len(x)], np.float64)
  yfit = np.zeros([len(x)], np.float64)

  for b in range(BootCyc):
    for i in range(len(x)):

      # Sample with/without replacement?
      if WithRep:
        j = np.random.randint(len(x))
      else:
        j = i

      # Sampling Statistical Uncertainty
      if not WithUncert or xsem[j] == 0.0:
        xtmp[i] = x[j]
      else:
        xtmp[i] = np.random.normal(x[j],xsem[j])
      if not WithUncert or ysem[j] == 0.0:
        ytmp[i] = y[j]
      else:
        ytmp[i] = np.random.normal(y[j],ysem[j])

      # Collect numbers for each boot cycle
      MBoot[0:Nm,b] = geterrormetrics(xtmp, ytmp)

  # Calculate mean and SEMs
  for m in range(Nm):
    MVals[m]=np.mean(MBoot[m])
    MSEMs[m]=np.std(MBoot[m])

  # Return Arrays
  return (MVals,MSEMs,MBoot)



### Load experimental file and place data in array
### Can't use np.loadtxt because there may be variable number of data points
### Print header thing
print "### Experimental file:",expfile,"ExpType:",exptype,"CorrectOA:",CorrectOA,"CorrectCB:",CorrectCB,"OAHOnly:",OAHOnly,"OAMeOnly:",OAMeOnly,"CalcPairDiffs:",CalcPairDiffs,"WithUncert:",WithUncert,"WithRep:",WithRep

with open(expfile, 'r') as expraw:
  explines = expraw.readlines()
exp = []
for line in explines:
  if not re.match(r'^\s*$', line):
    exp.append(line.rstrip().replace('\t', ' '))

### If flag, just do first six, or last six
if OAHOnly:
  exp = exp[0:6]
if OAMeOnly:
  exp = exp[6:12]
if OAHOnly and OAMeOnly:
  print "OAHOnly=True and OAMeOnly=True! Not compatible"
  exit()
N = len(exp)    # Number of data points in expfile
Np = (N-1)*N/2  # Number of data pairs

### Are these binding constants or free energy (or enthalpy)? Convert.
emean = np.zeros([N], np.float64)
esem = np.zeros([N], np.float64)
for i in range(len(exp)):
  cols = np.asarray(exp[i].split(), dtype=np.float64)
  if exptype == 'ka':          # Raw association constant
    dG = -R*T*np.log(cols)
    emean[i] = np.mean(dG)
    esem[i] = np.std(dG, ddof=1)/np.sqrt(len(dG))
  elif exptype == 'dg':        # Free energy at 298 K, kcal/mol, presumably
    emean[i] = cols[0]
    esem[i] = cols[1]
  elif exptype == 'dh':        # Enthalpy in calories, convert to kcal/mol
    emean[i] = np.mean(cols)/1000
    esem[i] = np.std(cols, ddof=1)/np.sqrt(len(cols))/1000
  else:
    print exptype, "... is not a valid experimental type"

### Calculate Experimental Pairwise Differences
if CalcPairDiffs:
  h = 0
  epmean = np.zeros([Np], np.float64)
  epsem = np.zeros([Np], np.float64)
  for i in range(len(exp)):
    for j in range(i+1, len(exp)):
      epmean[h] = emean[i] - emean[j]
      epsem[h] = np.sqrt( esem[i]**2 + esem[j]**2 )
      h += 1

### Setup arrays
### I should add a check to make sure number of data points is the same as experiment
Nc = len(calcfiles)
cmean = np.zeros([Nc,N], np.float64)
csem = np.zeros([Nc,N], np.float64)
if CalcPairDiffs:
  cdmean = np.zeros([Nc,Nd], np.float64)
  cdsem = np.zeros([Nc,Nd], np.float64)

RawMs = np.zeros([Nc,Nm], np.float64)
AllMBoot = np.zeros([Nc,Nm,BootCyc], np.float64)
AllMVals = np.zeros([Nc,Nm], np.float64)
AllMSEMs = np.zeros([Nc,Nm], np.float64)

### Print Header
#print "%30s " % ("Submission"),
print "%20s\t" % ("Submission"),
for Name in MNames:
  #print "%-22s " % (Name),
  print "%-s\tMean\t%1s\tSEM\t" % (Name,u'\u00B1'),
print ""

### Get/print numbers for each calculation method
nc = 0
for calcfile in calcfiles:
  calc = np.loadtxt(calcfile, np.float64)
  if OAHOnly:
    calc = calc[0:6]
  if OAMeOnly:
    calc = calc[6:12]
  for i in range(len(calc)):
    if np.isscalar(calc[i]) == True:      ### If scalar instead of array; ie, no SEM given.
      cmean[nc,i] = np.mean(calc[i])
      csem[nc,i] = 0.0
    else:                                 ### Assume Mean and SEM given
      cmean[nc,i] = calc[i,0]
      csem[nc,i] = calc[i,1]
    #print emean[i],esem[i],tmpcmean[i],tmpcsem[i]


  # Correct data set with MSE
  if CorrectOA:
    cmean[nc,0:6] = cmean[nc,0:6] - (np.mean(cmean[nc,0:6]) - np.mean(emean[0:6]))
    if len(calc) == 12:
      cmean[nc,6:12] = cmean[nc,6:12]- (np.mean(cmean[nc,6:12]) - np.mean(emean[6:12]))
  if CorrectCB:
    cmean[nc,0:10] = cmean[nc,0:10] - (np.mean(cmean[nc,0:10]) - np.mean(emean[0:10]))

  # Pairwise differences?
  if CalcPairDiffs:
    h = 0
    for i in range(len(calc)):
      for j in range(i+1, len(calc)):
        cpmean[nc,h] = cmean[nc,i] - cmean[nc,j]
        cpsem[nc,h] = np.sqrt( csem[nc,i]**2 + csem[nc,j]**2 )
        #print cpmean[nc,h]#,cdsem[nc,h]
        h += 1

  # Get raw error metrics (without resampling)
  RawMs[nc] = geterrormetrics(emean, cmean[nc])
  
  # Get resampled error metrics
  AllMVals[nc],AllMSEMs[nc],AllMBoot[nc] = bootstrap(emean, esem, cmean[nc], csem[nc])

  # Print Stuff
  #print "%30s " % (CalcNames[nc]),
  print "%20s\t" % (CalcNames[nc]),
  for m in range(Nm):
    #print "%6.2f(%6.2f %1s%6.2f) " % (RawMs[nc,m],AllMVals[nc,m],u'\u00B1',AllMSEMs[nc,m]),
    print "%6.2f\t%6.2f\t%1s\t%6.2f\t" % (RawMs[nc,m],AllMVals[nc,m],u'\u00B1',AllMSEMs[nc,m]),
  print ""
  sys.stdout.flush()

  nc += 1


