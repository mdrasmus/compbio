



text = os.popen("egrep '^generate:' ../spidir.log | sed 's/generate://'").read()
vals = map(float, text.split())

x = vals[::2]
y = vals[1::2]
#y2 = [exceptDefault(lambda i: exp(i), 0.0) for i in y]




p = plot(x, y, style='points')

#ests = [1.329375, 1.473125]
ests = map(float, os.popen("grep est_ ../spidir.log | awk '{print $2}'").read().split())

for est in ests[:5]:
    p.plot([est, est], [0, max(y)], style='lines')
