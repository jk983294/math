# colon operator, step is 1
print(5:8)  # 5 6 7 8
print(6.6:8.6)  # 6.6 7.6 8.6
# If the final element specified does not belong to the sequence then it is
# discarded.
print(3.8:8.4)  # 3.8 4.8 5.8 6.8 7.8

# seq operator, user defined step
seq(5, 9, by = 0.4)  # 5.0 5.4 5.8 6.2 6.6 7.0 7.4 7.8 8.2 8.6 9.0

# c() function, non-character values are coerced to character type if one of the
# elements is a character c stand for concatenation
v <- c("apple", "red", 5, TRUE)
print(v)  # !!! warn: different type convert to same type
