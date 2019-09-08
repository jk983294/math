# key is number
switch(1, "first", "second", "third")  # 'first'
switch(3, "first", "second", "third")  # 'third'
switch(0, "first", "second", "third")  # NULL
switch(4, "first", "second", "third")  # NULL

# key is string
switch("a", a = "x", b = "y", c = 5)  # 'x'
switch("c", a = "x", b = "y", c = 5)  # 5
switch("not_exist", a = "x", b = "y", c = 5)  # NULL

# another example
feelings <- c("sad", "afraid")
for (i in feelings) print(switch(i, happy = "I am glad you are happy", afraid = "There is nothing to fear", 
    sad = "Cheer up", angry = "Calm down now"))
