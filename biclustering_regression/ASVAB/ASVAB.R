
# Set working directory
# setwd()


new_data <- read.table('ASVAB.dat', sep=' ')
names(new_data) <- c('R0000100',
'R0536300',
'R0536401',
'R0536402',
'R1235800',
'R1318200',
'R1482600',
'R3961900',
'R3989200',
'R5473700',
'R7237400',
'R9705200',
'R9705300',
'R9705400',
'R9705500',
'R9705600',
'R9705700',
'R9705800',
'R9705900',
'R9706000',
'R9706100',
'R9706200',
'R9706300',
'R9706400',
'R9706500',
'R9706600',
'R9706700',
'R9706800',
'R9706900',
'R9707000',
'R9707100',
'R9707200',
'R9707300',
'R9707400',
'R9707500',
'R9793400',
'R9793500',
'R9793600',
'R9793700',
'R9793800',
'R9793900',
'S1552700')


# Handle missing values

new_data[new_data == -1] = NA  # Refused
new_data[new_data == -2] = NA  # Dont know
new_data[new_data == -3] = NA  # Invalid missing
new_data[new_data == -4] = NA  # Valid missing
new_data[new_data == -5] = NA  # Non-interview


# If there are values not categorized they will be represented as NA

vallabels = function(data) {
  data$R0536300 <- factor(data$R0536300,
levels=c(0.0,1.0,2.0),
labels=c("No Information",
"Male",
"Female"))
  data$R0536401 <- factor(data$R0536401,
levels=c(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0),
labels=c("1: January",
"2: February",
"3: March",
"4: April",
"5: May",
"6: June",
"7: July",
"8: August",
"9: September",
"10: October",
"11: November",
"12: December"))
  data$R1235800 <- factor(data$R1235800,
levels=c(0.0,1.0),
labels=c("Oversample",
"Cross-sectional"))
  data$R1482600 <- factor(data$R1482600,
levels=c(1.0,2.0,3.0,4.0),
labels=c("Black",
"Hispanic",
"Mixed Race (Non-Hispanic)",
"Non-Black / Non-Hispanic"))
return(data)
}


# If there are values not categorized they will be represented as NA

vallabels_continuous = function(data) {
data$R0000100[1.0 <= data$R0000100 & data$R0000100 <= 999.0] <- 1.0
data$R0000100[1000.0 <= data$R0000100 & data$R0000100 <= 1999.0] <- 1000.0
data$R0000100[2000.0 <= data$R0000100 & data$R0000100 <= 2999.0] <- 2000.0
data$R0000100[3000.0 <= data$R0000100 & data$R0000100 <= 3999.0] <- 3000.0
data$R0000100[4000.0 <= data$R0000100 & data$R0000100 <= 4999.0] <- 4000.0
data$R0000100[5000.0 <= data$R0000100 & data$R0000100 <= 5999.0] <- 5000.0
data$R0000100[6000.0 <= data$R0000100 & data$R0000100 <= 6999.0] <- 6000.0
data$R0000100[7000.0 <= data$R0000100 & data$R0000100 <= 7999.0] <- 7000.0
data$R0000100[8000.0 <= data$R0000100 & data$R0000100 <= 8999.0] <- 8000.0
data$R0000100[9000.0 <= data$R0000100 & data$R0000100 <= 9999.0] <- 9000.0
data$R0000100 <- factor(data$R0000100,
levels=c(0.0,1.0,1000.0,2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0),
labels=c("0",
"1 TO 999",
"1000 TO 1999",
"2000 TO 2999",
"3000 TO 3999",
"4000 TO 4999",
"5000 TO 5999",
"6000 TO 6999",
"7000 TO 7999",
"8000 TO 8999",
"9000 TO 9999"))
data$R1318200[50.0 <= data$R1318200 & data$R1318200 <= 59.0] <- 50.0
data$R1318200[60.0 <= data$R1318200 & data$R1318200 <= 69.0] <- 60.0
data$R1318200[70.0 <= data$R1318200 & data$R1318200 <= 79.0] <- 70.0
data$R1318200[80.0 <= data$R1318200 & data$R1318200 <= 89.0] <- 80.0
data$R1318200[90.0 <= data$R1318200 & data$R1318200 <= 99.0] <- 90.0
data$R1318200[100.0 <= data$R1318200 & data$R1318200 <= 109.0] <- 100.0
data$R1318200[110.0 <= data$R1318200 & data$R1318200 <= 119.0] <- 110.0
data$R1318200[120.0 <= data$R1318200 & data$R1318200 <= 129.0] <- 120.0
data$R1318200[130.0 <= data$R1318200 & data$R1318200 <= 139.0] <- 130.0
data$R1318200[140.0 <= data$R1318200 & data$R1318200 <= 149.0] <- 140.0
data$R1318200 <- factor(data$R1318200,
levels=c(50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0),
labels=c("50 TO 59",
"60 TO 69",
"70 TO 79",
"80 TO 89",
"90 TO 99",
"100 TO 109",
"110 TO 119",
"120 TO 129",
"130 TO 139",
"140 TO 149"))
data$R3961900[50.0 <= data$R3961900 & data$R3961900 <= 59.0] <- 50.0
data$R3961900[60.0 <= data$R3961900 & data$R3961900 <= 69.0] <- 60.0
data$R3961900[70.0 <= data$R3961900 & data$R3961900 <= 79.0] <- 70.0
data$R3961900[80.0 <= data$R3961900 & data$R3961900 <= 89.0] <- 80.0
data$R3961900[90.0 <= data$R3961900 & data$R3961900 <= 99.0] <- 90.0
data$R3961900[100.0 <= data$R3961900 & data$R3961900 <= 109.0] <- 100.0
data$R3961900[110.0 <= data$R3961900 & data$R3961900 <= 119.0] <- 110.0
data$R3961900[120.0 <= data$R3961900 & data$R3961900 <= 129.0] <- 120.0
data$R3961900[130.0 <= data$R3961900 & data$R3961900 <= 139.0] <- 130.0
data$R3961900[140.0 <= data$R3961900 & data$R3961900 <= 149.0] <- 140.0
data$R3961900 <- factor(data$R3961900,
levels=c(50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0),
labels=c("50 TO 59",
"60 TO 69",
"70 TO 79",
"80 TO 89",
"90 TO 99",
"100 TO 109",
"110 TO 119",
"120 TO 129",
"130 TO 139",
"140 TO 149"))
data$R3989200[50.0 <= data$R3989200 & data$R3989200 <= 59.0] <- 50.0
data$R3989200[60.0 <= data$R3989200 & data$R3989200 <= 69.0] <- 60.0
data$R3989200[70.0 <= data$R3989200 & data$R3989200 <= 79.0] <- 70.0
data$R3989200[80.0 <= data$R3989200 & data$R3989200 <= 89.0] <- 80.0
data$R3989200[90.0 <= data$R3989200 & data$R3989200 <= 99.0] <- 90.0
data$R3989200[100.0 <= data$R3989200 & data$R3989200 <= 109.0] <- 100.0
data$R3989200[110.0 <= data$R3989200 & data$R3989200 <= 119.0] <- 110.0
data$R3989200[120.0 <= data$R3989200 & data$R3989200 <= 129.0] <- 120.0
data$R3989200[130.0 <= data$R3989200 & data$R3989200 <= 139.0] <- 130.0
data$R3989200[140.0 <= data$R3989200 & data$R3989200 <= 149.0] <- 140.0
data$R3989200 <- factor(data$R3989200,
levels=c(50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0),
labels=c("50 TO 59",
"60 TO 69",
"70 TO 79",
"80 TO 89",
"90 TO 99",
"100 TO 109",
"110 TO 119",
"120 TO 129",
"130 TO 139",
"140 TO 149"))
data$R5473700[50.0 <= data$R5473700 & data$R5473700 <= 59.0] <- 50.0
data$R5473700[60.0 <= data$R5473700 & data$R5473700 <= 69.0] <- 60.0
data$R5473700[70.0 <= data$R5473700 & data$R5473700 <= 79.0] <- 70.0
data$R5473700[80.0 <= data$R5473700 & data$R5473700 <= 89.0] <- 80.0
data$R5473700[90.0 <= data$R5473700 & data$R5473700 <= 99.0] <- 90.0
data$R5473700[100.0 <= data$R5473700 & data$R5473700 <= 109.0] <- 100.0
data$R5473700[110.0 <= data$R5473700 & data$R5473700 <= 119.0] <- 110.0
data$R5473700[120.0 <= data$R5473700 & data$R5473700 <= 129.0] <- 120.0
data$R5473700[130.0 <= data$R5473700 & data$R5473700 <= 139.0] <- 130.0
data$R5473700[140.0 <= data$R5473700 & data$R5473700 <= 149.0] <- 140.0
data$R5473700 <- factor(data$R5473700,
levels=c(50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0),
labels=c("50 TO 59",
"60 TO 69",
"70 TO 79",
"80 TO 89",
"90 TO 99",
"100 TO 109",
"110 TO 119",
"120 TO 129",
"130 TO 139",
"140 TO 149"))
data$R7237400[50.0 <= data$R7237400 & data$R7237400 <= 59.0] <- 50.0
data$R7237400[60.0 <= data$R7237400 & data$R7237400 <= 69.0] <- 60.0
data$R7237400[70.0 <= data$R7237400 & data$R7237400 <= 79.0] <- 70.0
data$R7237400[80.0 <= data$R7237400 & data$R7237400 <= 89.0] <- 80.0
data$R7237400[90.0 <= data$R7237400 & data$R7237400 <= 99.0] <- 90.0
data$R7237400[100.0 <= data$R7237400 & data$R7237400 <= 109.0] <- 100.0
data$R7237400[110.0 <= data$R7237400 & data$R7237400 <= 119.0] <- 110.0
data$R7237400[120.0 <= data$R7237400 & data$R7237400 <= 129.0] <- 120.0
data$R7237400[130.0 <= data$R7237400 & data$R7237400 <= 139.0] <- 130.0
data$R7237400[140.0 <= data$R7237400 & data$R7237400 <= 149.0] <- 140.0
data$R7237400 <- factor(data$R7237400,
levels=c(50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0),
labels=c("50 TO 59",
"60 TO 69",
"70 TO 79",
"80 TO 89",
"90 TO 99",
"100 TO 109",
"110 TO 119",
"120 TO 129",
"130 TO 139",
"140 TO 149"))
data$R9705200[0.0 <= data$R9705200 & data$R9705200 <= 500.0] <- 0.0
data$R9705200[501.0 <= data$R9705200 & data$R9705200 <= 1000.0] <- 501.0
data$R9705200[1001.0 <= data$R9705200 & data$R9705200 <= 1500.0] <- 1001.0
data$R9705200[1501.0 <= data$R9705200 & data$R9705200 <= 2000.0] <- 1501.0
data$R9705200[2001.0 <= data$R9705200 & data$R9705200 <= 2500.0] <- 2001.0
data$R9705200[2501.0 <= data$R9705200 & data$R9705200 <= 3000.0] <- 2501.0
data$R9705200[3001.0 <= data$R9705200 & data$R9705200 <= 3500.0] <- 3001.0
data$R9705200 <- factor(data$R9705200,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9705300[0.0 <= data$R9705300 & data$R9705300 <= 500.0] <- 0.0
data$R9705300[501.0 <= data$R9705300 & data$R9705300 <= 1000.0] <- 501.0
data$R9705300[1001.0 <= data$R9705300 & data$R9705300 <= 1500.0] <- 1001.0
data$R9705300[1501.0 <= data$R9705300 & data$R9705300 <= 2000.0] <- 1501.0
data$R9705300[2001.0 <= data$R9705300 & data$R9705300 <= 2500.0] <- 2001.0
data$R9705300[2501.0 <= data$R9705300 & data$R9705300 <= 3000.0] <- 2501.0
data$R9705300[3001.0 <= data$R9705300 & data$R9705300 <= 3500.0] <- 3001.0
data$R9705300 <- factor(data$R9705300,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9705400[0.0 <= data$R9705400 & data$R9705400 <= 500.0] <- 0.0
data$R9705400[501.0 <= data$R9705400 & data$R9705400 <= 1000.0] <- 501.0
data$R9705400[1001.0 <= data$R9705400 & data$R9705400 <= 1500.0] <- 1001.0
data$R9705400[1501.0 <= data$R9705400 & data$R9705400 <= 2000.0] <- 1501.0
data$R9705400[2001.0 <= data$R9705400 & data$R9705400 <= 2500.0] <- 2001.0
data$R9705400[2501.0 <= data$R9705400 & data$R9705400 <= 3000.0] <- 2501.0
data$R9705400[3001.0 <= data$R9705400 & data$R9705400 <= 3500.0] <- 3001.0
data$R9705400 <- factor(data$R9705400,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9705500[0.0 <= data$R9705500 & data$R9705500 <= 500.0] <- 0.0
data$R9705500[501.0 <= data$R9705500 & data$R9705500 <= 1000.0] <- 501.0
data$R9705500[1001.0 <= data$R9705500 & data$R9705500 <= 1500.0] <- 1001.0
data$R9705500[1501.0 <= data$R9705500 & data$R9705500 <= 2000.0] <- 1501.0
data$R9705500[2001.0 <= data$R9705500 & data$R9705500 <= 2500.0] <- 2001.0
data$R9705500[2501.0 <= data$R9705500 & data$R9705500 <= 3000.0] <- 2501.0
data$R9705500[3001.0 <= data$R9705500 & data$R9705500 <= 3500.0] <- 3001.0
data$R9705500 <- factor(data$R9705500,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9705600[0.0 <= data$R9705600 & data$R9705600 <= 5000.0] <- 0.0
data$R9705600[5001.0 <= data$R9705600 & data$R9705600 <= 10000.0] <- 5001.0
data$R9705600[10001.0 <= data$R9705600 & data$R9705600 <= 15000.0] <- 10001.0
data$R9705600[15001.0 <= data$R9705600 & data$R9705600 <= 20000.0] <- 15001.0
data$R9705600[20001.0 <= data$R9705600 & data$R9705600 <= 25000.0] <- 20001.0
data$R9705600[25001.0 <= data$R9705600 & data$R9705600 <= 30000.0] <- 25001.0
data$R9705600[30001.0 <= data$R9705600 & data$R9705600 <= 35000.0] <- 30001.0
data$R9705600[35001.0 <= data$R9705600 & data$R9705600 <= 40000.0] <- 35001.0
data$R9705600[40001.0 <= data$R9705600 & data$R9705600 <= 45000.0] <- 40001.0
data$R9705600[45001.0 <= data$R9705600 & data$R9705600 <= 50000.0] <- 45001.0
data$R9705600 <- factor(data$R9705600,
levels=c(0.0,5001.0,10001.0,15001.0,20001.0,25001.0,30001.0,35001.0,40001.0,45001.0),
labels=c("0 TO 5000",
"5001 TO 10000",
"10001 TO 15000",
"15001 TO 20000",
"20001 TO 25000",
"25001 TO 30000",
"30001 TO 35000",
"35001 TO 40000",
"40001 TO 45000",
"45001 TO 50000"))
data$R9705700[0.0 <= data$R9705700 & data$R9705700 <= 5000.0] <- 0.0
data$R9705700[5001.0 <= data$R9705700 & data$R9705700 <= 10000.0] <- 5001.0
data$R9705700[10001.0 <= data$R9705700 & data$R9705700 <= 15000.0] <- 10001.0
data$R9705700[15001.0 <= data$R9705700 & data$R9705700 <= 20000.0] <- 15001.0
data$R9705700[20001.0 <= data$R9705700 & data$R9705700 <= 25000.0] <- 20001.0
data$R9705700[25001.0 <= data$R9705700 & data$R9705700 <= 30000.0] <- 25001.0
data$R9705700[30001.0 <= data$R9705700 & data$R9705700 <= 35000.0] <- 30001.0
data$R9705700[35001.0 <= data$R9705700 & data$R9705700 <= 40000.0] <- 35001.0
data$R9705700[40001.0 <= data$R9705700 & data$R9705700 <= 45000.0] <- 40001.0
data$R9705700[45001.0 <= data$R9705700 & data$R9705700 <= 50000.0] <- 45001.0
data$R9705700 <- factor(data$R9705700,
levels=c(0.0,5001.0,10001.0,15001.0,20001.0,25001.0,30001.0,35001.0,40001.0,45001.0),
labels=c("0 TO 5000",
"5001 TO 10000",
"10001 TO 15000",
"15001 TO 20000",
"20001 TO 25000",
"25001 TO 30000",
"30001 TO 35000",
"35001 TO 40000",
"40001 TO 45000",
"45001 TO 50000"))
data$R9705800[0.0 <= data$R9705800 & data$R9705800 <= 500.0] <- 0.0
data$R9705800[501.0 <= data$R9705800 & data$R9705800 <= 1000.0] <- 501.0
data$R9705800[1001.0 <= data$R9705800 & data$R9705800 <= 1500.0] <- 1001.0
data$R9705800[1501.0 <= data$R9705800 & data$R9705800 <= 2000.0] <- 1501.0
data$R9705800[2001.0 <= data$R9705800 & data$R9705800 <= 2500.0] <- 2001.0
data$R9705800[2501.0 <= data$R9705800 & data$R9705800 <= 3000.0] <- 2501.0
data$R9705800[3001.0 <= data$R9705800 & data$R9705800 <= 3500.0] <- 3001.0
data$R9705800 <- factor(data$R9705800,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9705900[0.0 <= data$R9705900 & data$R9705900 <= 500.0] <- 0.0
data$R9705900[501.0 <= data$R9705900 & data$R9705900 <= 1000.0] <- 501.0
data$R9705900[1001.0 <= data$R9705900 & data$R9705900 <= 1500.0] <- 1001.0
data$R9705900[1501.0 <= data$R9705900 & data$R9705900 <= 2000.0] <- 1501.0
data$R9705900[2001.0 <= data$R9705900 & data$R9705900 <= 2500.0] <- 2001.0
data$R9705900[2501.0 <= data$R9705900 & data$R9705900 <= 3000.0] <- 2501.0
data$R9705900[3001.0 <= data$R9705900 & data$R9705900 <= 3500.0] <- 3001.0
data$R9705900 <- factor(data$R9705900,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9706000[0.0 <= data$R9706000 & data$R9706000 <= 500.0] <- 0.0
data$R9706000[501.0 <= data$R9706000 & data$R9706000 <= 1000.0] <- 501.0
data$R9706000[1001.0 <= data$R9706000 & data$R9706000 <= 1500.0] <- 1001.0
data$R9706000[1501.0 <= data$R9706000 & data$R9706000 <= 2000.0] <- 1501.0
data$R9706000[2001.0 <= data$R9706000 & data$R9706000 <= 2500.0] <- 2001.0
data$R9706000[2501.0 <= data$R9706000 & data$R9706000 <= 3000.0] <- 2501.0
data$R9706000[3001.0 <= data$R9706000 & data$R9706000 <= 3500.0] <- 3001.0
data$R9706000 <- factor(data$R9706000,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9706100[0.0 <= data$R9706100 & data$R9706100 <= 500.0] <- 0.0
data$R9706100[501.0 <= data$R9706100 & data$R9706100 <= 1000.0] <- 501.0
data$R9706100[1001.0 <= data$R9706100 & data$R9706100 <= 1500.0] <- 1001.0
data$R9706100[1501.0 <= data$R9706100 & data$R9706100 <= 2000.0] <- 1501.0
data$R9706100[2001.0 <= data$R9706100 & data$R9706100 <= 2500.0] <- 2001.0
data$R9706100[2501.0 <= data$R9706100 & data$R9706100 <= 3000.0] <- 2501.0
data$R9706100[3001.0 <= data$R9706100 & data$R9706100 <= 3500.0] <- 3001.0
data$R9706100 <- factor(data$R9706100,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9706200[0.0 <= data$R9706200 & data$R9706200 <= 500.0] <- 0.0
data$R9706200[501.0 <= data$R9706200 & data$R9706200 <= 1000.0] <- 501.0
data$R9706200[1001.0 <= data$R9706200 & data$R9706200 <= 1500.0] <- 1001.0
data$R9706200[1501.0 <= data$R9706200 & data$R9706200 <= 2000.0] <- 1501.0
data$R9706200[2001.0 <= data$R9706200 & data$R9706200 <= 2500.0] <- 2001.0
data$R9706200[2501.0 <= data$R9706200 & data$R9706200 <= 3000.0] <- 2501.0
data$R9706200[3001.0 <= data$R9706200 & data$R9706200 <= 3500.0] <- 3001.0
data$R9706200 <- factor(data$R9706200,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9706300[0.0 <= data$R9706300 & data$R9706300 <= 500.0] <- 0.0
data$R9706300[501.0 <= data$R9706300 & data$R9706300 <= 1000.0] <- 501.0
data$R9706300[1001.0 <= data$R9706300 & data$R9706300 <= 1500.0] <- 1001.0
data$R9706300[1501.0 <= data$R9706300 & data$R9706300 <= 2000.0] <- 1501.0
data$R9706300[2001.0 <= data$R9706300 & data$R9706300 <= 2500.0] <- 2001.0
data$R9706300[2501.0 <= data$R9706300 & data$R9706300 <= 3000.0] <- 2501.0
data$R9706300[3001.0 <= data$R9706300 & data$R9706300 <= 3500.0] <- 3001.0
data$R9706300 <- factor(data$R9706300,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9706400[0.0 <= data$R9706400 & data$R9706400 <= 500.0] <- 0.0
data$R9706400[501.0 <= data$R9706400 & data$R9706400 <= 1000.0] <- 501.0
data$R9706400[1001.0 <= data$R9706400 & data$R9706400 <= 1500.0] <- 1001.0
data$R9706400[1501.0 <= data$R9706400 & data$R9706400 <= 2000.0] <- 1501.0
data$R9706400[2001.0 <= data$R9706400 & data$R9706400 <= 2500.0] <- 2001.0
data$R9706400[2501.0 <= data$R9706400 & data$R9706400 <= 3000.0] <- 2501.0
data$R9706400[3001.0 <= data$R9706400 & data$R9706400 <= 3500.0] <- 3001.0
data$R9706400 <- factor(data$R9706400,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9706500[0.0 <= data$R9706500 & data$R9706500 <= 500.0] <- 0.0
data$R9706500[501.0 <= data$R9706500 & data$R9706500 <= 1000.0] <- 501.0
data$R9706500[1001.0 <= data$R9706500 & data$R9706500 <= 1500.0] <- 1001.0
data$R9706500[1501.0 <= data$R9706500 & data$R9706500 <= 2000.0] <- 1501.0
data$R9706500[2001.0 <= data$R9706500 & data$R9706500 <= 2500.0] <- 2001.0
data$R9706500[2501.0 <= data$R9706500 & data$R9706500 <= 3000.0] <- 2501.0
data$R9706500[3001.0 <= data$R9706500 & data$R9706500 <= 3500.0] <- 3001.0
data$R9706500 <- factor(data$R9706500,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9706600[0.0 <= data$R9706600 & data$R9706600 <= 500.0] <- 0.0
data$R9706600[501.0 <= data$R9706600 & data$R9706600 <= 1000.0] <- 501.0
data$R9706600[1001.0 <= data$R9706600 & data$R9706600 <= 1500.0] <- 1001.0
data$R9706600[1501.0 <= data$R9706600 & data$R9706600 <= 2000.0] <- 1501.0
data$R9706600[2001.0 <= data$R9706600 & data$R9706600 <= 2500.0] <- 2001.0
data$R9706600[2501.0 <= data$R9706600 & data$R9706600 <= 3000.0] <- 2501.0
data$R9706600[3001.0 <= data$R9706600 & data$R9706600 <= 3500.0] <- 3001.0
data$R9706600 <- factor(data$R9706600,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9706700[0.0 <= data$R9706700 & data$R9706700 <= 500.0] <- 0.0
data$R9706700[501.0 <= data$R9706700 & data$R9706700 <= 1000.0] <- 501.0
data$R9706700[1001.0 <= data$R9706700 & data$R9706700 <= 1500.0] <- 1001.0
data$R9706700[1501.0 <= data$R9706700 & data$R9706700 <= 2000.0] <- 1501.0
data$R9706700[2001.0 <= data$R9706700 & data$R9706700 <= 2500.0] <- 2001.0
data$R9706700[2501.0 <= data$R9706700 & data$R9706700 <= 3000.0] <- 2501.0
data$R9706700[3001.0 <= data$R9706700 & data$R9706700 <= 3500.0] <- 3001.0
data$R9706700 <- factor(data$R9706700,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9706800[0.0 <= data$R9706800 & data$R9706800 <= 1000.0] <- 0.0
data$R9706800[1001.0 <= data$R9706800 & data$R9706800 <= 2000.0] <- 1001.0
data$R9706800[2001.0 <= data$R9706800 & data$R9706800 <= 3000.0] <- 2001.0
data$R9706800[3001.0 <= data$R9706800 & data$R9706800 <= 4000.0] <- 3001.0
data$R9706800[4001.0 <= data$R9706800 & data$R9706800 <= 5000.0] <- 4001.0
data$R9706800[5001.0 <= data$R9706800 & data$R9706800 <= 6000.0] <- 5001.0
data$R9706800[6001.0 <= data$R9706800 & data$R9706800 <= 7000.0] <- 6001.0
data$R9706800[7001.0 <= data$R9706800 & data$R9706800 <= 8000.0] <- 7001.0
data$R9706800[8001.0 <= data$R9706800 & data$R9706800 <= 9000.0] <- 8001.0
data$R9706800[9001.0 <= data$R9706800 & data$R9706800 <= 10000.0] <- 9001.0
data$R9706800 <- factor(data$R9706800,
levels=c(0.0,1001.0,2001.0,3001.0,4001.0,5001.0,6001.0,7001.0,8001.0,9001.0),
labels=c("0 TO 1000",
"1001 TO 2000",
"2001 TO 3000",
"3001 TO 4000",
"4001 TO 5000",
"5001 TO 6000",
"6001 TO 7000",
"7001 TO 8000",
"8001 TO 9000",
"9001 TO 10000"))
data$R9706900[0.0 <= data$R9706900 & data$R9706900 <= 1000.0] <- 0.0
data$R9706900[1001.0 <= data$R9706900 & data$R9706900 <= 2000.0] <- 1001.0
data$R9706900[2001.0 <= data$R9706900 & data$R9706900 <= 3000.0] <- 2001.0
data$R9706900[3001.0 <= data$R9706900 & data$R9706900 <= 4000.0] <- 3001.0
data$R9706900[4001.0 <= data$R9706900 & data$R9706900 <= 5000.0] <- 4001.0
data$R9706900[5001.0 <= data$R9706900 & data$R9706900 <= 6000.0] <- 5001.0
data$R9706900[6001.0 <= data$R9706900 & data$R9706900 <= 7000.0] <- 6001.0
data$R9706900[7001.0 <= data$R9706900 & data$R9706900 <= 8000.0] <- 7001.0
data$R9706900[8001.0 <= data$R9706900 & data$R9706900 <= 9000.0] <- 8001.0
data$R9706900[9001.0 <= data$R9706900 & data$R9706900 <= 10000.0] <- 9001.0
data$R9706900 <- factor(data$R9706900,
levels=c(0.0,1001.0,2001.0,3001.0,4001.0,5001.0,6001.0,7001.0,8001.0,9001.0),
labels=c("0 TO 1000",
"1001 TO 2000",
"2001 TO 3000",
"3001 TO 4000",
"4001 TO 5000",
"5001 TO 6000",
"6001 TO 7000",
"7001 TO 8000",
"8001 TO 9000",
"9001 TO 10000"))
data$R9707000[0.0 <= data$R9707000 & data$R9707000 <= 500.0] <- 0.0
data$R9707000[501.0 <= data$R9707000 & data$R9707000 <= 1000.0] <- 501.0
data$R9707000[1001.0 <= data$R9707000 & data$R9707000 <= 1500.0] <- 1001.0
data$R9707000[1501.0 <= data$R9707000 & data$R9707000 <= 2000.0] <- 1501.0
data$R9707000[2001.0 <= data$R9707000 & data$R9707000 <= 2500.0] <- 2001.0
data$R9707000[2501.0 <= data$R9707000 & data$R9707000 <= 3000.0] <- 2501.0
data$R9707000[3001.0 <= data$R9707000 & data$R9707000 <= 3500.0] <- 3001.0
data$R9707000 <- factor(data$R9707000,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9707100[0.0 <= data$R9707100 & data$R9707100 <= 500.0] <- 0.0
data$R9707100[501.0 <= data$R9707100 & data$R9707100 <= 1000.0] <- 501.0
data$R9707100[1001.0 <= data$R9707100 & data$R9707100 <= 1500.0] <- 1001.0
data$R9707100[1501.0 <= data$R9707100 & data$R9707100 <= 2000.0] <- 1501.0
data$R9707100[2001.0 <= data$R9707100 & data$R9707100 <= 2500.0] <- 2001.0
data$R9707100[2501.0 <= data$R9707100 & data$R9707100 <= 3000.0] <- 2501.0
data$R9707100[3001.0 <= data$R9707100 & data$R9707100 <= 3500.0] <- 3001.0
data$R9707100 <- factor(data$R9707100,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9707200[0.0 <= data$R9707200 & data$R9707200 <= 500.0] <- 0.0
data$R9707200[501.0 <= data$R9707200 & data$R9707200 <= 1000.0] <- 501.0
data$R9707200[1001.0 <= data$R9707200 & data$R9707200 <= 1500.0] <- 1001.0
data$R9707200[1501.0 <= data$R9707200 & data$R9707200 <= 2000.0] <- 1501.0
data$R9707200[2001.0 <= data$R9707200 & data$R9707200 <= 2500.0] <- 2001.0
data$R9707200[2501.0 <= data$R9707200 & data$R9707200 <= 3000.0] <- 2501.0
data$R9707200[3001.0 <= data$R9707200 & data$R9707200 <= 3500.0] <- 3001.0
data$R9707200 <- factor(data$R9707200,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9707300[0.0 <= data$R9707300 & data$R9707300 <= 500.0] <- 0.0
data$R9707300[501.0 <= data$R9707300 & data$R9707300 <= 1000.0] <- 501.0
data$R9707300[1001.0 <= data$R9707300 & data$R9707300 <= 1500.0] <- 1001.0
data$R9707300[1501.0 <= data$R9707300 & data$R9707300 <= 2000.0] <- 1501.0
data$R9707300[2001.0 <= data$R9707300 & data$R9707300 <= 2500.0] <- 2001.0
data$R9707300[2501.0 <= data$R9707300 & data$R9707300 <= 3000.0] <- 2501.0
data$R9707300[3001.0 <= data$R9707300 & data$R9707300 <= 3500.0] <- 3001.0
data$R9707300 <- factor(data$R9707300,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9707400[0.0 <= data$R9707400 & data$R9707400 <= 500.0] <- 0.0
data$R9707400[501.0 <= data$R9707400 & data$R9707400 <= 1000.0] <- 501.0
data$R9707400[1001.0 <= data$R9707400 & data$R9707400 <= 1500.0] <- 1001.0
data$R9707400[1501.0 <= data$R9707400 & data$R9707400 <= 2000.0] <- 1501.0
data$R9707400[2001.0 <= data$R9707400 & data$R9707400 <= 2500.0] <- 2001.0
data$R9707400[2501.0 <= data$R9707400 & data$R9707400 <= 3000.0] <- 2501.0
data$R9707400[3001.0 <= data$R9707400 & data$R9707400 <= 3500.0] <- 3001.0
data$R9707400 <- factor(data$R9707400,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9707500[0.0 <= data$R9707500 & data$R9707500 <= 500.0] <- 0.0
data$R9707500[501.0 <= data$R9707500 & data$R9707500 <= 1000.0] <- 501.0
data$R9707500[1001.0 <= data$R9707500 & data$R9707500 <= 1500.0] <- 1001.0
data$R9707500[1501.0 <= data$R9707500 & data$R9707500 <= 2000.0] <- 1501.0
data$R9707500[2001.0 <= data$R9707500 & data$R9707500 <= 2500.0] <- 2001.0
data$R9707500[2501.0 <= data$R9707500 & data$R9707500 <= 3000.0] <- 2501.0
data$R9707500[3001.0 <= data$R9707500 & data$R9707500 <= 3500.0] <- 3001.0
data$R9707500 <- factor(data$R9707500,
levels=c(0.0,501.0,1001.0,1501.0,2001.0,2501.0,3001.0),
labels=c("0 TO 500",
"501 TO 1000",
"1001 TO 1500",
"1501 TO 2000",
"2001 TO 2500",
"2501 TO 3000",
"3001 TO 3500"))
data$R9793400[0.0 <= data$R9793400 & data$R9793400 <= 10.0] <- 0.0
data$R9793400[11.0 <= data$R9793400 & data$R9793400 <= 20.0] <- 11.0
data$R9793400[21.0 <= data$R9793400 & data$R9793400 <= 30.0] <- 21.0
data$R9793400[31.0 <= data$R9793400 & data$R9793400 <= 40.0] <- 31.0
data$R9793400[41.0 <= data$R9793400 & data$R9793400 <= 50.0] <- 41.0
data$R9793400 <- factor(data$R9793400,
levels=c(0.0,11.0,21.0,31.0,41.0),
labels=c("0 TO 10",
"11 TO 20",
"21 TO 30",
"31 TO 40",
"41 TO 50"))
data$R9793500[0.0 <= data$R9793500 & data$R9793500 <= 10.0] <- 0.0
data$R9793500[11.0 <= data$R9793500 & data$R9793500 <= 20.0] <- 11.0
data$R9793500[21.0 <= data$R9793500 & data$R9793500 <= 30.0] <- 21.0
data$R9793500[31.0 <= data$R9793500 & data$R9793500 <= 40.0] <- 31.0
data$R9793500[41.0 <= data$R9793500 & data$R9793500 <= 50.0] <- 41.0
data$R9793500 <- factor(data$R9793500,
levels=c(0.0,11.0,21.0,31.0,41.0),
labels=c("0 TO 10",
"11 TO 20",
"21 TO 30",
"31 TO 40",
"41 TO 50"))
data$R9793600[0.0 <= data$R9793600 & data$R9793600 <= 10.0] <- 0.0
data$R9793600[11.0 <= data$R9793600 & data$R9793600 <= 20.0] <- 11.0
data$R9793600[21.0 <= data$R9793600 & data$R9793600 <= 30.0] <- 21.0
data$R9793600[31.0 <= data$R9793600 & data$R9793600 <= 40.0] <- 31.0
data$R9793600[41.0 <= data$R9793600 & data$R9793600 <= 50.0] <- 41.0
data$R9793600 <- factor(data$R9793600,
levels=c(0.0,11.0,21.0,31.0,41.0),
labels=c("0 TO 10",
"11 TO 20",
"21 TO 30",
"31 TO 40",
"41 TO 50"))
data$R9793700[0.0 <= data$R9793700 & data$R9793700 <= 10.0] <- 0.0
data$R9793700[11.0 <= data$R9793700 & data$R9793700 <= 20.0] <- 11.0
data$R9793700[21.0 <= data$R9793700 & data$R9793700 <= 30.0] <- 21.0
data$R9793700[31.0 <= data$R9793700 & data$R9793700 <= 40.0] <- 31.0
data$R9793700[41.0 <= data$R9793700 & data$R9793700 <= 50.0] <- 41.0
data$R9793700 <- factor(data$R9793700,
levels=c(0.0,11.0,21.0,31.0,41.0),
labels=c("0 TO 10",
"11 TO 20",
"21 TO 30",
"31 TO 40",
"41 TO 50"))
data$R9793800[1.0 <= data$R9793800 & data$R9793800 <= 100.0] <- 1.0
data$R9793800[101.0 <= data$R9793800 & data$R9793800 <= 200.0] <- 101.0
data$R9793800[201.0 <= data$R9793800 & data$R9793800 <= 300.0] <- 201.0
data$R9793800[301.0 <= data$R9793800 & data$R9793800 <= 400.0] <- 301.0
data$R9793800[401.0 <= data$R9793800 & data$R9793800 <= 500.0] <- 401.0
data$R9793800[501.0 <= data$R9793800 & data$R9793800 <= 600.0] <- 501.0
data$R9793800[601.0 <= data$R9793800 & data$R9793800 <= 700.0] <- 601.0
data$R9793800[701.0 <= data$R9793800 & data$R9793800 <= 800.0] <- 701.0
data$R9793800 <- factor(data$R9793800,
levels=c(1.0,101.0,201.0,301.0,401.0,501.0,601.0,701.0),
labels=c("1 TO 100",
"101 TO 200",
"201 TO 300",
"301 TO 400",
"401 TO 500",
"501 TO 600",
"601 TO 700",
"701 TO 800"))
data$R9793900[1.0 <= data$R9793900 & data$R9793900 <= 100.0] <- 1.0
data$R9793900[101.0 <= data$R9793900 & data$R9793900 <= 200.0] <- 101.0
data$R9793900[201.0 <= data$R9793900 & data$R9793900 <= 300.0] <- 201.0
data$R9793900[301.0 <= data$R9793900 & data$R9793900 <= 400.0] <- 301.0
data$R9793900[401.0 <= data$R9793900 & data$R9793900 <= 500.0] <- 401.0
data$R9793900[501.0 <= data$R9793900 & data$R9793900 <= 600.0] <- 501.0
data$R9793900[601.0 <= data$R9793900 & data$R9793900 <= 700.0] <- 601.0
data$R9793900[701.0 <= data$R9793900 & data$R9793900 <= 800.0] <- 701.0
data$R9793900 <- factor(data$R9793900,
levels=c(1.0,101.0,201.0,301.0,401.0,501.0,601.0,701.0),
labels=c("1 TO 100",
"101 TO 200",
"201 TO 300",
"301 TO 400",
"401 TO 500",
"501 TO 600",
"601 TO 700",
"701 TO 800"))
data$S1552700[50.0 <= data$S1552700 & data$S1552700 <= 59.0] <- 50.0
data$S1552700[60.0 <= data$S1552700 & data$S1552700 <= 69.0] <- 60.0
data$S1552700[70.0 <= data$S1552700 & data$S1552700 <= 79.0] <- 70.0
data$S1552700[80.0 <= data$S1552700 & data$S1552700 <= 89.0] <- 80.0
data$S1552700[90.0 <= data$S1552700 & data$S1552700 <= 99.0] <- 90.0
data$S1552700[100.0 <= data$S1552700 & data$S1552700 <= 109.0] <- 100.0
data$S1552700[110.0 <= data$S1552700 & data$S1552700 <= 119.0] <- 110.0
data$S1552700[120.0 <= data$S1552700 & data$S1552700 <= 129.0] <- 120.0
data$S1552700[130.0 <= data$S1552700 & data$S1552700 <= 139.0] <- 130.0
data$S1552700[140.0 <= data$S1552700 & data$S1552700 <= 149.0] <- 140.0
data$S1552700 <- factor(data$S1552700,
levels=c(50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0),
labels=c("50 TO 59",
"60 TO 69",
"70 TO 79",
"80 TO 89",
"90 TO 99",
"100 TO 109",
"110 TO 119",
"120 TO 129",
"130 TO 139",
"140 TO 149"))
return(data)
}

varlabels <- c("PUBID - YTH ID CODE 1997",
"KEY!SEX (SYMBOL) 1997",
"KEY!BDATE M/Y (SYMBOL) 1997",
"KEY!BDATE M/Y (SYMBOL) 1997",
"CV_SAMPLE_TYPE 1997",
"CV_PIAT_STANDARD_UPD 1997",
"KEY!RACE_ETHNICITY (SYMBOL) 1997",
"CV_PIAT_STANDARD_UPD 1999",
"CV_PIAT_STANDARD_UPD 1998",
"CV_PIAT_STANDARD_SCORE 2000",
"CV_PIAT_STANDARD_SCORE 2001",
"ASVAB_GS_ABILITY_EST_POS 1999",
"ASVAB_AR_ABILITY_EST_POS 1999",
"ASVAB_WK_ABILITY_EST_POS 1999",
"ASVAB_PC_ABILITY_EST_POS 1999",
"ASVAB_NO_ABILITY_EST_POS 1999",
"ASVAB_CS_ABILITY_EST_POS 1999",
"ASVAB_AI_ABILITY_EST_POS 1999",
"ASVAB_SI_ABILITY_EST_POS 1999",
"ASVAB_MK_ABILITY_EST_POS 1999",
"ASVAB_MC_ABILITY_EST_POS 1999",
"ASVAB_EI_ABILITY_EST_POS 1999",
"ASVAB_AO_ABILITY_EST_POS 1999",
"ASVAB_GS_ABILITY_EST_NEG 1999",
"ASVAB_AR_ABILITY_EST_NEG 1999",
"ASVAB_WK_ABILITY_EST_NEG 1999",
"ASVAB_PC_ABILITY_EST_NEG 1999",
"ASVAB_NO_ABILITY_EST_NEG 1999",
"ASVAB_CS_ABILITY_EST_NEG 1999",
"ASVAB_AI_ABILITY_EST_NEG 1999",
"ASVAB_SI_ABILITY_EST_NEG 1999",
"ASVAB_MK_ABILITY_EST_NEG 1999",
"ASVAB_MC_ABILITY_EST_NEG 1999",
"ASVAB_EI_ABILITY_EST_NEG 1999",
"ASVAB_AO_ABILITY_EST_NEG 1999",
"TRANS_ACT_COMP HSTR",
"TRANS_ACT_ENG HSTR",
"TRANS_ACT_MATH HSTR",
"TRANS_ACT_READ HSTR",
"TRANS_SAT_VERBAL HSTR",
"TRANS_SAT_MATH HSTR",
"CV_PIAT_STANDARD_SCORE 2002"
)


# Use qnames rather than rnums

qnames = function(data) {
names(data) <- c("PUBID_1997",
"KEY_SEX_1997",
"KEY_BDATE_M_1997",
"KEY_BDATE_Y_1997",
"CV_SAMPLE_TYPE_1997",
"CV_PIAT_STANDARD_UPD_1997",
"KEY_RACE_ETHNICITY_1997",
"CV_PIAT_STANDARD_UPD_1999",
"CV_PIAT_STANDARD_UPD_1998",
"CV_PIAT_STANDARD_SCORE_2000",
"CV_PIAT_STANDARD_SCORE_2001",
"ASVAB_GS_ABILITY_EST_POS_1999",
"ASVAB_AR_ABILITY_EST_POS_1999",
"ASVAB_WK_ABILITY_EST_POS_1999",
"ASVAB_PC_ABILITY_EST_POS_1999",
"ASVAB_NO_ABILITY_EST_POS_1999",
"ASVAB_CS_ABILITY_EST_POS_1999",
"ASVAB_AI_ABILITY_EST_POS_1999",
"ASVAB_SI_ABILITY_EST_POS_1999",
"ASVAB_MK_ABILITY_EST_POS_1999",
"ASVAB_MC_ABILITY_EST_POS_1999",
"ASVAB_EI_ABILITY_EST_POS_1999",
"ASVAB_AO_ABILITY_EST_POS_1999",
"ASVAB_GS_ABILITY_EST_NEG_1999",
"ASVAB_AR_ABILITY_EST_NEG_1999",
"ASVAB_WK_ABILITY_EST_NEG_1999",
"ASVAB_PC_ABILITY_EST_NEG_1999",
"ASVAB_NO_ABILITY_EST_NEG_1999",
"ASVAB_CS_ABILITY_EST_NEG_1999",
"ASVAB_AI_ABILITY_EST_NEG_1999",
"ASVAB_SI_ABILITY_EST_NEG_1999",
"ASVAB_MK_ABILITY_EST_NEG_1999",
"ASVAB_MC_ABILITY_EST_NEG_1999",
"ASVAB_EI_ABILITY_EST_NEG_1999",
"ASVAB_AO_ABILITY_EST_NEG_1999",
"TRANS_ACT_COMP_HSTR",
"TRANS_ACT_ENG_HSTR",
"TRANS_ACT_MATH_HSTR",
"TRANS_ACT_READ_HSTR",
"TRANS_SAT_VERBAL_HSTR",
"TRANS_SAT_MATH_HSTR",
"CV_PIAT_STANDARD_SCORE_2002")
return(data)
}


#********************************************************************************************************

# Remove the '#' before the following line to create a data file called "categories" with value labels.
#categories <- vallabels(new_data)

# Remove the '#' before the following lines to rename variables using Qnames instead of Reference Numbers
#new_data <- qnames(new_data)
#categories <- qnames(categories)

# Produce summaries for the raw (uncategorized) data file
summary(new_data)

# Remove the '#' before the following lines to produce summaries for the "categories" data file.
#categories <- vallabels(new_data)
#categories <- vallabels_continuous(new_data)
#summary(categories)

#************************************************************************************************************

