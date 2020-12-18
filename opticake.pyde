# area of polygon with vertices
def area(poly):
    return 0.5 * abs(sum(poly[i][0] * poly[i - 1][1] - poly[i - 1][0] * poly[i][1] for i in xrange(len(poly))))

def line_points(a, b, n):
    return [(lerp(a[0], b[0], float(i) / n), lerp(a[1], b[1], float(i) / n)) for i in xrange(n + 1)]

theta0 = 60
x0 = sqrt(1 / (4 * (1 + tan(radians(theta0 - 45)))))
y0 = x0 + x0 * tan(radians(theta0 - 45))
N = 6
mutdev = 0.0002
# get area correct to a lot of decimal places
area_penalty = 10 ** 10
coarse_tol = 0.03
sN = 100
total_mutations = 0

print(x0,  y0)

l1 = line_points((y0, 0), (x0, x0), N)

def draw_line(l):
    for i in xrange(len(l) - 1):
        ellipse(width * l[i][0], height * l[i][1], 3, 3)
        line(width * l[i][0], height * l[i][1], width * l[i + 1][0], height * l[i + 1][1])
    ellipse(width * l[i + 1][0], height * l[i + 1][1], 3, 3)

def line_length(l):
    return sum(dist(l[i][0], l[i][1], l[i + 1][0], l[i + 1][1]) for i in xrange(len(l) - 1))

def objective_fn(l1, l2, l3, l4, l5):
    # total_length = line_length(l1) + line_length(l1[-1:] + l2) + line_length(l3) + line_length(l3[-1:] + l4) + line_length(l1[-1:] + l5 + l3[-1:])
    total_length = 4 * line_length(l1) + dist(l1[-1][0], l1[-1][1], 1 - l1[-1][0], 1 - l1[-1][1])
    A1 = area(l1 + l2 + [(0, 0)])
    # A2 = area(l3 + l4 + [(1, 1)])
    # A3 = area(l3 + l5[::-1] + l1[-1:] + l2 + [(0, 1)])
    A4 = area(l1 + l5 + l3[-1:] + l4 + [(1, 0)])
    return (total_length + 2 * area_penalty * (abs(A1 - 0.25) ** 2 + abs(A4 - 0.25) ** 2),
        total_length, A1, A1, A4, A4)

# for theta0 in [50 + i * 20 / float(sN) for i in range(sN)]:
#     x0 = sqrt(1 / (4 * (1 + tan(radians(theta0 - 45)))))
#     y0 = x0 + x0 * tan(radians(theta0 - 45))
#     tl1 = line_points((y0, 0), (x0, x0), N)
#     tl2 = [(y, x) for x, y in tl1][-1::-1]
#     tl3 = [(1 - x, 1 - y) for x, y in tl1]
#     tl4 = [(y, x) for x, y in tl3][-1::-1]
#     tl5 = []
#     print(objective_fn(tl1, tl2, tl3, tl4, tl5))

def mutate_coord(x):
    if x == 0:
        return x
    return x + randomGaussian() * mutdev

def mutate_line(l):
    return [(mutate_coord(x), mutate_coord(y)) if ind != len(l) - 1 else (mutate_coord(x),) * 2 for ind, (x, y) in enumerate(l)]

def mutate(cur_obval, l1, l2, l3, l4, l5):
    global total_mutations
    cand_l1 = mutate_line(l1)
    # cand_l2 = mutate_line(l2)
    # cand_l3 = mutate_line(l3)
    # cand_l4 = mutate_line(l4)
    # cand_l5 = mutate_line(l5)
    cand_l2 = [(y, x) for x, y in cand_l1][-1::-1]
    cand_l3 = [(1 - x, 1 - y) for x, y in cand_l1]
    cand_l4 = [(y, x) for x, y in cand_l3][-1::-1]
    cand_l5 = []
    obval, _, _, _, _, _ = objective_fn(cand_l1, cand_l2, cand_l3, cand_l4, cand_l5)
    if obval < cur_obval:
        l1[:] = cand_l1
        # l2[:] = cand_l2
        # l3[:] = cand_l3
        # l4[:] = cand_l4
        # l5[:] = cand_l5
        total_mutations += 1
        return obval
    return cur_obval

def recoarse(l):
    skip_next = False
    for i in xrange(len(l) - 3, -1, -1):
        if skip_next:
            skip_next = False
            continue
        if dist(l[i][0], l[i][1], l[i + 1][0], l[i + 1][1]) < coarse_tol:
            skip_next = True
            del l[i + 1]

def refine(l):
    new_l = []
    for i in xrange(len(l) - 1):
        new_l.append(l[i])
        new_l.append((0.5 * (l[i][0] + l[i + 1][0]), 0.5 * (l[i][1] + l[i + 1][1])))
    new_l.append(l[-1])
    l[:] = new_l

def setup():
    global f
    size(800, 800)
    ellipseMode(RADIUS)
    f = createFont("courier", 14)

def draw():
    background(0)
    stroke(255)
    textFont(f, 14)
    l2 = [(y, x) for x, y in l1][-1::-1]
    l3 = [(1 - x, 1 - y) for x, y in l1]
    l4 = [(y, x) for x, y in l3][-1::-1]
    l5 = []
    draw_line(l1)
    draw_line(l1[-1:] + l2)
    draw_line(l3)
    draw_line(l3[-1:] + l4)
    draw_line(l1[-1:] + l5 + l3[-1:])
    obval, tlen, A1, A2, A3, A4 = objective_fn(l1, l2, l3, l4, l5)
    text("len: {}".format(tlen), 10, 24)
    text("obj: {}".format(obval), 10, 40)
    text(" A1: {}".format(A1), 10, 56)
    text(" A2: {}".format(A2), 10, 72)
    text(" A3: {}".format(A3), 10, 88)
    text(" A4: {}".format(A4), 10, 104)
    text("fps: {}".format(frameRate), 10, 120)
    text("mut: {}".format(mutdev), 10, 136)
    text("pnt: {0} ({0:.0e})".format(area_penalty), 10, 152)
    text("tot: {}".format(total_mutations), 10, 168)
    for _ in xrange(500):
        obval = mutate(obval, l1, l2, l3, l4, l5)

def keyPressed():
    global mutdev, area_penalty
    if keyCode == ord("C"):
        print("coarse")
        recoarse(l1)
    elif keyCode == ord("F"):
        print("fine")
        refine(l1)
    elif keyCode == ord("P"):
        l2 = [(y, x) for x, y in l1][-1::-1]
        l3 = [(1 - x, 1 - y) for x, y in l1]
        l4 = [(y, x) for x, y in l3][-1::-1]
        l5 = []
        print(objective_fn(l1, l2, l3, l4, l5))
        print("area penalty: {}".format(area_penalty))
        print("{} points".format(len(l1)))
        print(l1)
    elif keyCode == ord("M"):
        if key == "m":
            mutdev /= 2
        if key == "M":
            mutdev *= 2
    elif keyCode == ord("A"):
        if key == "a":
            area_penalty /= 10
        if key == "A":
            area_penalty *= 10
