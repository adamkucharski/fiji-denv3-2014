# Load data non-zero weekly reports
# Github version

timeser = timeser14

colID=(names(timeser)==locationtab[iiH])

endwk=length(timeser[,1])
startdate=min(as.Date(timeser$date,origin="1970-01-01"))
enddate=max(as.Date(timeser$date,origin="1970-01-01"))

pickdate=(as.Date(timeser$date,origin="1970-01-01")>=startdate & as.Date(timeser$date,origin="1970-01-01")<=enddate)
y.vals = as.numeric(timeser[,colID])[pickdate]

# Choose first data point above zero
firstentry=min(c(1:endwk)[(y.vals>0)])
pickdate=c(firstentry:endwk)

y.vals=y.vals[firstentry:endwk]

# Add data onto end ** Need to start data sets at 11-04 for seasonality **
date_list=as.Date(timeser$date,origin="1970-01-01")[pickdate]
time.vals = as.numeric(date_list-min(date_list)) #+ shift_to_nov  # Shift by start date so both time series line up -- IMPORTANT FOR SEASONALITY
time.vals = time.vals+time.vals[2]-time.vals[1] # adjust if first value non zero

shift_to_nov = (min(date_list) - start.date - 7) %>% as.numeric() # Adjust so seasonality matches #(firstentry-1) * 7 # No longer needed as start date fixed

add.null.dates = ((sample.collection.date - ( min(date_list) + max(time.vals) ) )/7) %>% ceiling() %>% as.numeric()  # Up to start of serosurvey -- 16/10/15
#add.null.dates = ((as.Date("2014-09-22") - ( min(date_list) + max(time.vals) ) )/7) %>% ceiling() %>% as.numeric()  # Up to start of serosurvey -- 16/10/15

#add.null.dates = 10 # Only fit to end point

time.vals = c(time.vals,seq(max(time.vals)+7,max(time.vals)+7*add.null.dates,7)) # Add extra onto end if needed

rR.vals = 1 # Set reporting to one
rep.drop.date = sum(date_list<=swap.date) #Set up date of drop

