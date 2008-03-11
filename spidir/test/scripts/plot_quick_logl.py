



logls = map2(float, readDelim("logl.txt"))


if 1:
    plot(cget(logls, 1), cget(logls, 2))

if 0:
    logls.sort(key=lambda x: x[2], reverse=True)
    plot(cget(logls, 0))

if 0:
    logls.sort(key=lambda x: x[0], reverse=True)
    plot(cget(logls, 2))

