# Binary Files Read the 'mtcars' data frame as a csv file and store only the
# columns 'cyl','am' and 'gear'
write.table(mtcars, file = "mtcars.csv", row.names = FALSE, na = "", col.names = TRUE, 
    sep = ",")

# Store 5 records from the csv file as a new data frame
new.mtcars <- read.table("mtcars.csv", sep = ",", header = TRUE, nrows = 5)

# Create a connection object to write the binary file using mode 'wb'
write.filename = file("/web/com/binmtcars.dat", "wb")

# Write the column names of the data frame to the connection object
writeBin(colnames(new.mtcars), write.filename)

# Write the records in each of the column to the file
writeBin(c(new.mtcars$cyl, new.mtcars$am, new.mtcars$gear), write.filename)

# Close the file for writing so that it can be read by other program
close(write.filename)

# Create a connection object to read the file in binary mode using 'rb'
read.filename <- file("/web/com/binmtcars.dat", "rb")

# First read the column names. n=3 as we have 3 columns
column.names <- readBin(read.filename, character(), n = 3)

# Next read the column values. n=18 as we have 3 column names and 15 values
read.filename <- file("/web/com/binmtcars.dat", "rb")
bindata <- readBin(read.filename, integer(), n = 18)

# Read the values from 4th byte to 8th byte which represents 'cyl'
cyldata = bindata[4:8]

# Read the values form 9th byte to 13th byte which represents 'am'
amdata = bindata[9:13]

# Read the values form 9th byte to 13th byte which represents 'gear'
geardata = bindata[14:18]

# Combine all the read values to a dat frame
finaldata = cbind(cyldata, amdata, geardata)
colnames(finaldata) = column.names
print(finaldata)
