

x = []
vals = []

for line in os.popen("grep sample_int ../tmp/out.log | sed 's/sample_int: //'"):
    tokens = line.split()
    x.append(int(tokens[0]))
    vals.append(float(tokens[1]))

s=1100
t=800
plot(x[s:s+t], vals[s:s+t], style='lines')



