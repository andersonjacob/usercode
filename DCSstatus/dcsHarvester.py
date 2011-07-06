#!/usr/bin/python

import cx_Oracle
import datetime
import os

from xml.dom import minidom
from time import gmtime, strftime
from optparse import OptionParser

import sys
#import string

scriptLocation = sys.argv[0]
scriptDir = '/'
scriptDir = scriptDir.join(scriptLocation.split('/')[:-1])

if scriptLocation == scriptDir:
    scriptDir = '.'

if len(scriptDir) < 1:
    scriptDir = '.'
#print scriptLocation,scriptDir

tns_admin = ''
selftns = scriptDir
if os.environ.has_key('TNS_ADMIN'):
    tns_admin = os.environ['TNS_ADMIN']
else:
    tns_admin = None

if os.access(selftns, os.R_OK):
    os.environ['TNS_ADMIN'] = selftns
else:
    os.environ['TNS_ADMIN'] = '/nfshome0/andersj/oracle'

# hardwired lumisection length
OneLS = datetime.timedelta(seconds=23.31041)

#parse command line options
parser = OptionParser()
parser.add_option("-i", "--init", action="store_true", dest="init",
                  default=False, help="action to create init file")
parser.add_option("-c", "--current", action="store_true", dest="current",
                  default=False, help="use current values in the last value"+ \
                  " table to create the init file, mutually exclusive with "+ \
                  "--init")
parser.add_option("-u", "--update", action="store_true", dest="update",
                  default=False, help="action to create update file")
parser.add_option("-r", "--run", type="int", dest="run",
                  help="required : run number for actions", metavar="RUN")
parser.add_option("--keep", action="store_true", dest="keep", default=False,
                  help="keep xml don't zip; implies --nocopy")
parser.add_option("--nocopy", action="store_true", dest="nocopy",
                  default=False, help="don't copy zip file to spool dir")
parser.add_option("--offset", type="int", dest="offset", 
                  default=0, metavar="OFFSET", 
                  help="offset for the sub-version; normally 0")
#parser.add_option("--nodev", action="store_false", dest="dev", default=True,
#                  help="use OMDS instead of the develpment DB")
parser.add_option("--dev", action="store_true", dest="dev", default=False,
                  help="spool to the development DB instead of OMDS")
parser.add_option("--directory", dest="outd",
                  help="directory where the xml and zip files will be " +\
                  "created and stored before spooling.  " +\
                  "Defaults to " + scriptDir, metavar="DIR",
                  default=scriptDir)
parser.add_option("--tty", dest="tty", action="store_true", default=False,
                  help="use tty output")

(options, args) = parser.parse_args()

if (options.run == None) or (options.run < 1):
    parser.error("invalid run number: a run number is required")

if (options.init) and (options.current):
    parser.error("init and current options are mutually exclusive.")

if options.keep:
    options.nocopy = True

logfilename = '/tmp/dcsHarvester.log'
if (options.tty):
    logfilename = '/dev/tty'

sys.stdout = open(logfilename, 'a', 1)
sys.stderr = open(logfilename, 'a', 1)

print "TNS_ADMIN:",os.environ['TNS_ADMIN']
print "file createion path:",options.outd
if not os.access(options.outd, os.F_OK):
    os.makedirs(options.outd)
#print options
#print args

# connect to the database
connectString = "CMS_HCL_APPUSER_R/HCAL_Reader_44@"

## if (options.dev):
##     connectString += "cmsdevr_lb"
## else:
connectString += "cms_omds_lb"

## #this is kluging some things together to use the dev db
## if (options.dev) and (tns_admin != '/nfshome0/andersj/oracle'):
##     print "Correcting TNS_ADMIN..."
##     os.environ['TNS_ADMIN'] = '/nfshome0/andersj/oracle'
##     print "TNS_ADMIN:", os.environ['TNS_ADMIN']

omds = cx_Oracle.connect(connectString)
cursor = omds.cursor()

#RunNumber = 123596

# function to create a new text node for the xml
def newTextNode(xml, nodeName, nodeText):
    newNode = xml.createElement(nodeName)
    newNode.appendChild(xml.createTextNode(str(nodeText)))
    return newNode

# function to make a python datetime from the oracle Timestamp
def ToDatetime(dt):
    return datetime.datetime(dt.year, dt.month, dt.day, dt.hour, dt.minute, \
                             dt.second, dt.fsecond)

# determine the length in seconds of a timedelta
def secLen(dt) :
    return dt.days*24*3600 + dt.seconds + dt.microseconds*0.000001

# determine the lumi section number given a timedelta and the lumi section 
# length as another timedelta
def LSNum(dt, LSlen):
    return int(secLen(dt)//secLen(LSlen) + 2)

#create a header node for the XML file
def createXMLHeader(xml, runNum, initUpdate):
    header = xml.createElement("HEADER")
    htype = xml.createElement("TYPE")
    htype.appendChild(newTextNode(xml, "EXTENSION_TABLE_NAME",
                                  "HCAL_DCS_ENV_VALUES_V1"))
    htype.appendChild(newTextNode(xml, "NAME",
                                  "Hcal DCS Env Values [" + initUpdate + "]"))
    run = xml.createElement("RUN")
    run.appendChild(newTextNode(xml, "RUN_TYPE", "DCS Voltage Values"))
    run.appendChild(newTextNode(xml, "RUN_NUMBER", runNum))
    header.appendChild(htype)
    header.appendChild(run)
    return header

def createElements(xml, runNum, dev=True):
    elements = xml.createElement("ELEMENTS")
    dset = xml.createElement("DATA_SET")
    dset.setAttribute("id", "-1")
    elements.appendChild(dset)
    iov = xml.createElement("IOV")
    iov.setAttribute("id", "1")
    iov.appendChild(newTextNode(xml, "INTERVAL_OF_VALIDITY_BEGIN", runNum))
    iov.appendChild(newTextNode(xml, "INTERVAL_OF_VALIDITY_END", -1))
    elements.appendChild(iov)
    tag = xml.createElement("TAG")
    tag.setAttribute("id", "2")
    tag.setAttribute("mode", "auto")
    tagName = "HcalDcsValues_v1.00_"
    if (dev):
        tagName += "test"
    else:
        tagName += "offline"
    tag.appendChild(newTextNode(xml, "TAG_NAME", tagName))
    tag.appendChild(newTextNode(xml, "DETECTOR_NAME", "HCAL"))
    tag.appendChild(newTextNode(xml, "COMMENT_DESCRIPTION",
                                "Run " + str(runNum) + " DCS values"))
    elements.appendChild(iov)
    elements.appendChild(tag)
    return elements

def createMaps(xml):
    maps = xml.createElement("MAPS")
    tag = xml.createElement("TAG")
    tag.setAttribute("idref", "2")
    iov = xml.createElement("IOV")
    iov.setAttribute("idref", "1")
    dset = xml.createElement("DATA_SET")
    dset.setAttribute("idref", "-1")
    iov.appendChild(dset)
    tag.appendChild(iov)
    maps.appendChild(tag)
    return maps

#add the stuff to the dataset node that comes before the data
def addDataSetFrontMatter(xml, dataSet, runNum, tols, update):
    dataSet.appendChild(newTextNode(xml, "VERSION", "V_" + str(runNum) + \
                                        "_T_" + tols))
    dataSet.appendChild(newTextNode(xml, "SUBVERSION", update + options.offset))
    dataSet.appendChild(newTextNode(xml, "CREATED_BY_USER", "Jake Anderson"))
    dataSet.appendChild(newTextNode(xml, "CREATE_TIMESTAMP",
                                    strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    dataSet.appendChild(newTextNode(xml, "COMMENT_DESCRIPTION",
                                    "Voltages for run " + str(runNum)))

# use the status to determine the validity of the value and zero it out
# if needed.
def valueForStatus(value, status):
    if (status == 192):
        return value
    elif (status == -1):
        return 0.
    elif (status == 0):
        return 0.
    else :
        return value
    

#create a data node
def createDataNode(xml, dpname, value, changeTime, status, LS = -10):
    dataNode = xml.createElement("DATA")
    dataNode.appendChild(newTextNode(xml, "VALUE_TIMESTAMP", changeTime))
    dataNode.appendChild(newTextNode(xml, "VALUE",
                                     valueForStatus(value,status)))
    dataNode.appendChild(newTextNode(xml, "DPNAME", dpname))
    dataNode.appendChild(newTextNode(xml, "VALUE_STATUS", status))
    if (LS >= 0):
        dataNode.appendChild(newTextNode(xml, "LUMISECTION_BEGIN", LS))
    return dataNode

#add the standalone attribute to the xml declaration
def addStandalone(xml):
    xmlstring = xml.toprettyxml( indent="   ", encoding="UTF-8" , newl="\n")
    return xmlstring.replace('?>', ' standalone="yes"?>', 1)

#write updates by trying to catch communication glitches
def writeUpdate(xml, dataSet, channelData, T0):
    i = 0
    currV = 0.
    currS = 0
    lastS = 0
    while i < len(channelData):
        row = channelData[i]
        diff = row[3] - T0
        lsnum = LSNum(diff, OneLS)
        if lsnum < 1:
            lsnum = 1
        if (row[2] != None):
            currV = row[2]
        if (row[4] != None):
            currS = row[4]      
        #print row,lsnum,currV,currS,lastS
        if (currS != -1):
            dataSet.appendChild(createDataNode(xml, row[1], currV,
                                               row[3], currS, lsnum))
            lastS = currS
        else:
            nextS = -999
            if (i+1 < len(channelData)):
                nextrow = channelData[i+1]
                nextS = nextrow[4]
            if (nextS != lastS):
                dataSet.appendChild(createDataNode(xml, row[1], currV, row[3],
                                                   currS, lsnum))
                lastS = currS
            else:
                i+=1
        i+=1
    
#create the xml ouput file
def outputFile(T0, tols, update):

    xml = minidom.Document()
    root = xml.createElement("ROOT")
    root.setAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
    xml.appendChild(root)

    initUpdate = "Update"
    if (update == 0):
        initUpdate = "Init"

    root.appendChild(createXMLHeader(xml, options.run, initUpdate))

    dataSet = xml.createElement("DATA_SET")
    addDataSetFrontMatter(xml, dataSet, options.run, tols, update)
    
    root.appendChild(dataSet)

    currDPID = 0
    channelData = []

    for row in cursor.fetchall():
        if (update == 0):
            dataSet.appendChild(createDataNode(xml, row[1], row[2], T0, row[4]))
        else:
            if (currDPID > 0) and (currDPID != row[0]):
                writeUpdate(xml, dataSet, channelData, T0)
                currDPID = row[0]
                channelData = []
            channelData.append(row)
                                
            ## dpchange = ToDatetime(row[3])
            ## diff = dpchange - T0
            ## lsnum = LSNum(diff,OneLS)
            ## if lsnum < 1:
            ##     lsnum = 1
            ## dataSet.appendChild(createDataNode(xml,row[1],row[2],row[3],lsnum))
    if (update != 0):
        writeUpdate(xml, dataSet, channelData, T0)

        
    root.appendChild(createElements(xml, options.run, options.dev))
    root.appendChild(createMaps(xml))
    outFilename = None
    print
    if (cursor.rowcount > 0):
        outFilename = options.outd + "/Hcal_DCS_Values_Run" + \
                      str(options.run) + "_" + initUpdate + ".xml"
        outf = open(outFilename, "w")
        outf.write(addStandalone(xml))
        outf.close()

        if not options.keep:
            if (os.system("zip -DjmT " + \
                          outFilename.replace(".xml", ".zip") + \
                          " " + outFilename + " >> " + logfilename + " 2>&1") == 0):
                outFilename = outFilename.replace(".xml", ".zip")
            else:
                raise "error zipping file"
            
    return outFilename

#copy the xml file to the spool directory and then delete it if successful
def spoolOutput(fname):
    outDir = " /var/spool/xmlloader/hcal/"
    if options.dev:
        outDir += "dev/"
    else:
        outDir += "prod/"
    outDir += "conditions/"
    if (os.system("cp -v " + fname + outDir + " >> " + logfilename + " 2>&1") == 0):
        os.system("rm -v " + fname + " >> " + logfilename + " 2>&1")
    else:
        print "failed to spool " + fname + " to" + outDir +".  Are you " + \
              "using a machine with the xml loader?"

#execute the query and then call the file creation and spooling.
def createAndSpoolFile(sql, runStart, runEnd, tols, update, LS0end):
    params = {}
    if ((options.init) or (options.update)):
        params = {'runStart': runStart}
    if (update==1):
        params['runEnd'] = runEnd

    try:
        cursor.execute(sql, params)
        if (update == 0):
            outfilename = outputFile(T0=runStart, update=update, tols=tols)
        else:
            outfilename = outputFile(T0=LS0end, update=update, tols=tols)

        if outfilename:
            print outfilename + " file created."
            if not options.nocopy:
                spoolOutput(outfilename)
        else:
            print "no update records; no file created"
    except cx_Oracle.DatabaseError, exc:
        err, = exc.args
        print "Oracle error:", err.code, "\nError message:", err.message, \
              "\nError context:", err.context
        print cursor.statement
        raise "Database Error"
    
# query to get run information for a given run number
runsql = """
SELECT RSUM.RUNNUMBER, RSUM.STARTTIME, RSUM.STOPTIME
FROM CMS_WBM.RUNSUMMARY RSUM
WHERE RSUM.RUNNUMBER = :runNumber
"""

# get the run start and end information from the database
try:
    cursor.execute(runsql, runNumber=options.run)
except cx_Oracle.DatabaseError, exc:
    err, = exc.args
    print "Oracle error:", err.code, "\nError message:", err.message, \
          "\nError context:", err.context
    print cursor.statement
    raise "Database Error"

RunRow = cursor.fetchone()
print RunRow
#RunStart = ToDatetime(RunRow[1])
RunStart = RunRow[1]
try:
    #RunEnd = ToDatetime(RunRow[2])
    RunEnd = RunRow[2]
    if (RunEnd < RunStart):
        RunEnd = datetime.datetime.utcnow()
except AttributeError:
    # if there is a problem because the run stoptime is None we'll just use
    # the current time
    RunEnd = datetime.datetime.utcnow()
print "At", str(datetime.datetime.utcnow()) + ", processing run", \
      str(RunRow[0]), "that started on", str(RunStart), "..."

lslensql = """
select cms_gt.FUNC_GET_LS_SEC_BY_RUNNR(:runNumber)
from dual
"""

lssql = """
SELECT stop_time, lsnr
FROM CMS_GT_MON.GT_MON_LS_STOP_TIME_VIEW
WHERE (runnr = :runNumber)
  AND (lsnr < 2)
order by lsnr
"""

EndLS0 = RunStart+OneLS

if options.update:
    try:
        cursor.execute(lslensql, runNumber=options.run)
    except cx_Oracle.DatabaseError, exc:
        err, = exc.args
        print "Oracle error:", err.code, "\nError message:", err.message, \
              "\nError context:", err.context
        print cursor.statement
        raise "Database Error"

    RunRow = cursor.fetchone()
    OneLS = datetime.timedelta(seconds=RunRow[0])

    try:
        cursor.execute(lssql, runNumber=options.run)
    except cx_Oracle.DatabaseError, exc:
        err, = exc.args
        print "Oracle error:", err.code, "\nError message:", err.message, \
              "\nError context:", err.context
        print cursor.statement
        raise "Database Error"

    RunRow = cursor.fetchone()
    #EndLS0 = ToDatetime(RunRow[0])
    EndLS0 = RunRow[0]

print "End of first LS:",str(EndLS0),"LS len:",str(secLen(OneLS)),"seconds"

tolversql = """
select distinct cds.version, cds.record_insertion_time, tol.start_timestamp
from cms_hcl_core_cond.cond_data_sets cds, 
  cms_hcl_core_cond.kinds_of_conditions koc, 
  cms_hcl_hcal_cond.hcal_dcs_env_tolerances_v1 tol
where (cds.kind_of_condition_id = koc.kind_of_condition_id)
  and (cds.condition_data_set_id = tol.condition_data_set_id)
  and (cds.is_record_deleted = 'F')
  and (koc.is_record_deleted = 'F')
  and (koc.name = 'Hcal DCS Env Tolerances [type 1]')
  and (tol.start_timestamp <= :runStart)
order by tol.start_timestamp desc, cds.record_insertion_time desc
"""

tolVer = ""
try:
    cursor.execute(tolversql, runStart=RunStart)
except cx_Oracle.DatabaseError, exc:
    err, = exc.args
    print "Oracle error:", err.code, "\nError message:", err.message, \
        "\nError context:", err.context
    print cursor.statement
    raise "Database Error"

RunRow = cursor.fetchone()
tolVer = RunRow[0]


# query to get the current values
currentsql = """
SELECT HVCHAN.DPID, CHAN_NAME.DPNAME, HVCHAN.VALUE_NUMBER, HVCHAN.CHANGE_DATE, HVCHAN.CHAN_STATUS
FROM 
  (SELECT DPID,
    MAX(DECODE (DPE_NAME, 'ACTUAL_VREAD', VALUE_NUMBER)) VALUE_NUMBER,
    MAX(DECODE (DPE_NAME, 'ACTUAL_VREAD', CHANGE_DATE)) CHANGE_DATE,
    MAX(DECODE (DPE_NAME, 'ACTUAL_STATUS', VALUE_NUMBER)) CHAN_STATUS
    FROM
	CMS_HCL_CORE_PVSS_COND.HCALHVCHANNEL_LV HCALHV
    GROUP BY DPID) HVCHAN, CMS_HCL_CORE_PVSS_COND.DP_NAME2ID CHAN_NAME
WHERE (HVCHAN.DPID = CHAN_NAME.ID)
ORDER BY CHAN_NAME.DPNAME
"""

## replaced this query Feb 2011. --Jake
## currentsql = """
## SELECT DP.ID, DP.DPNAME, LAST_V.VALUE_NUMBER, LAST_V.CHANGE_DATE, NULL
## FROM CMS_HCL_CORE_PVSS_COND.HCALHVCHANNEL_LV LAST_V,
##   CMS_HCL_CORE_PVSS_COND.DP_NAME2ID DP
## WHERE (LAST_V.DPID = DP.ID)
##   AND (LAST_V.DPE_NAME = 'ACTUAL_VREAD')
## ORDER BY DP.DPNAME
## """

## Not going to actually do anything anymore --Jake Feb 2011
## if options.current:
##     createAndSpoolFile(sql=currentsql, runStart=RunStart, runEnd=RunEnd, 
##                        tols=tolVer, update=0, LS0end=EndLS0)

# query to get initial values at the start of a run specified by the
# run number and RunStart
initsql = """
SELECT DCSHV.ID DPID, DCSHV.DPNAME, DCSHV.ACTUAL_VREAD, DCSHV.CHANGE_DATE,
  DCSSTAT.ACTUAL_STATUS, DCSSTAT.CHANGE_DATE
FROM 
  (
  SELECT DP.ID, DP.DPNAME, HVCHAN.ACTUAL_VREAD, HVCHAN.CHANGE_DATE,
    HVCHAN.ACTUAL_STATUS,
    ROW_NUMBER() OVER (PARTITION BY ID ORDER BY UPDATEID DESC) RN
  FROM CMS_HCL_CORE_PVSS_COND.HCALHVCHANNEL HVCHAN,
    CMS_HCL_CORE_PVSS_COND.DP_NAME2ID DP
  WHERE (HVCHAN.DPID = DP.ID)
    AND (HVCHAN.ACTUAL_VREAD is not null)
    AND (HVCHAN.CHANGE_DATE <= :runStart)
  ) DCSHV,
  (
  SELECT DP.ID, DP.DPNAME, HVCHAN.ACTUAL_VREAD, HVCHAN.CHANGE_DATE,
    HVCHAN.ACTUAL_STATUS,
    ROW_NUMBER() OVER (PARTITION BY ID ORDER BY UPDATEID DESC) RN
  FROM CMS_HCL_CORE_PVSS_COND.HCALHVCHANNEL HVCHAN,
    CMS_HCL_CORE_PVSS_COND.DP_NAME2ID DP
  WHERE (HVCHAN.DPID = DP.ID)
    AND (HVCHAN.ACTUAL_STATUS is not null)
    AND (HVCHAN.CHANGE_DATE <= :runStart)
  ) DCSSTAT
WHERE (DCSHV.RN = 1)
  AND (DCSSTAT.RN = 1)
  AND (DCSHV.ID = DCSSTAT.ID)
ORDER BY DCSHV.DPNAME
"""

## old query replaced Feb 2011. --Jake
## initsql = """
## SELECT DCSHV.ID DPID, DCSHV.DPNAME, DCSHV.ACTUAL_VREAD, DCSHV.CHANGE_DATE,
##   DCSHV.ACTUAL_STATUS
## FROM (
##   SELECT DP.ID, DP.DPNAME, HVCHAN.ACTUAL_VREAD, HVCHAN.CHANGE_DATE,
##     HVCHAN.ACTUAL_STATUS,
##     ROW_NUMBER() OVER (PARTITION BY ID ORDER BY UPDATEID DESC) RN
##   FROM CMS_HCL_CORE_PVSS_COND.HCALHVCHANNEL HVCHAN,
##     CMS_HCL_CORE_PVSS_COND.DP_NAME2ID DP
##   WHERE (HVCHAN.DPID = DP.ID)
##     AND (HVCHAN.ACTUAL_VREAD is NOT NULL)
##     AND (HVCHAN.ACTUAL_STATUS is NULL)
##     AND (HVCHAN.CHANGE_DATE <= :runStart)
## ) DCSHV
## WHERE (DCSHV.RN = 1)
## ORDER BY DCSHV.DPNAME
## """

## Not going to actually do anything anymore --Jake Feb 2011
## # get the initial values for a run and create the output file
## if options.init:
##     createAndSpoolFile(sql=initsql, runStart=RunStart, runEnd=RunEnd, 
##                        tols=tolVer, update=0, LS0end=EndLS0)

# query to get the changes during a run
updatesSql = """
SELECT DP.ID DPID, DP.DPNAME, 
       DCSHV.ACTUAL_VREAD VREAD, DCSHV.CHANGE_DATE V_CHANGE_DATE,
       DCSSTAT.ACTUAL_STATUS STATUS, DCSSTAT.CHANGE_DATE S_CHANGE_DATE
FROM 
  (
  SELECT HVCHAN.DPID, HVCHAN.ACTUAL_VREAD, HVCHAN.CHANGE_DATE,
    ROW_NUMBER() OVER (PARTITION BY DPID ORDER BY UPDATEID DESC) RN
  FROM CMS_HCL_CORE_PVSS_COND.HCALHVCHANNEL HVCHAN
  WHERE (HVCHAN.ACTUAL_VREAD is not null)
    AND (HVCHAN.CHANGE_DATE <= :runStart)
  ) DCSHV,
  (
  SELECT HVCHAN.DPID, HVCHAN.CHANGE_DATE, HVCHAN.ACTUAL_STATUS,
    ROW_NUMBER() OVER (PARTITION BY DPID ORDER BY UPDATEID DESC) RN
  FROM CMS_HCL_CORE_PVSS_COND.HCALHVCHANNEL HVCHAN
  WHERE (HVCHAN.ACTUAL_STATUS is not null)
    AND (HVCHAN.CHANGE_DATE <= :runStart)
  ) DCSSTAT, CMS_HCL_CORE_PVSS_COND.DP_NAME2ID DP
WHERE (DCSHV.RN = 1)
  AND (DCSSTAT.RN = 1)
  AND (DCSHV.DPID = DCSSTAT.DPID)
  AND (DP.ID = DCSHV.DPID)
--ORDER BY DP.DPNAME, DCSHV.CHANGE_DATE
UNION
SELECT DP.ID DPID, DP.DPNAME, 
       HVCHAN.ACTUAL_VREAD, HVCHAN.CHANGE_DATE V_CHANGE_DATE, 
       HVCHAN.ACTUAL_STATUS, HVCHAN.CHANGE_DATE S_CHANGE_DATE
FROM CMS_HCL_CORE_PVSS_COND.HCALHVCHANNEL HVCHAN,
  CMS_HCL_CORE_PVSS_COND.DP_NAME2ID DP
WHERE (HVCHAN.DPID = DP.ID)
  AND (
    (HVCHAN.ACTUAL_VREAD is not null) or 
    (HVCHAN.ACTUAL_STATUS is not null)
  )
  AND (HVCHAN.CHANGE_DATE < :runEnd)
  AND (HVCHAN.CHANGE_DATE > :runStart)
ORDER BY DPNAME, V_CHANGE_DATE, S_CHANGE_DATE
"""

## old query --Jake
## updatesSql = """
## SELECT DP.ID DPID, DP.DPNAME, HVCHAN.ACTUAL_VREAD, HVCHAN.CHANGE_DATE, 
##   HVCHAN.ACTUAL_STATUS
## FROM CMS_HCL_CORE_PVSS_COND.HCALHVCHANNEL HVCHAN,
##   CMS_HCL_CORE_PVSS_COND.DP_NAME2ID DP
## WHERE (HVCHAN.DPID = DP.ID)
##   AND (HVCHAN.CHANGE_DATE < :runEnd)
##   AND (HVCHAN.CHANGE_DATE > :runStart)
## ORDER BY DP.DPNAME, HVCHAN.CHANGE_DATE
## """

# get the update values for a run and create the output file
if options.update:
    createAndSpoolFile(sql=updatesSql, runStart=RunStart, runEnd=RunEnd, 
                       tols=tolVer, update=1, LS0end=EndLS0)

if tns_admin == None:
    del os.environ['TNS_ADMIN']
else:
    os.environ['TNS_ADMIN'] = tns_admin

print
