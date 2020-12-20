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
# l1 = [(0.5218420474908216, 0.0), (0.5218275327447004, 0.00897247052471134), (0.5217586831047495, 0.017950433258809753), (0.5216463353225416, 0.026864199857383735), (0.5214915590646519, 0.03575230400097134), (0.5212913195959596, 0.04465102916553483), (0.5210399056106595, 0.053535962725398754), (0.5207492575814839, 0.062438139116591546), (0.5204114869703174, 0.07133637195781008), (0.5200355287063358, 0.08020274047079425), (0.5196007800833624, 0.08936649243071747), (0.5191259015691034, 0.09820854221492796), (0.5186054133329656, 0.1071376027638102), (0.5180420825484222, 0.11597323979383105), (0.5174401624107501, 0.12477871254507551), (0.5168028613156669, 0.1335154790810992), (0.5161173102555479, 0.14229467117051317), (0.5152629712823763, 0.1523360206837114), (0.5143574717240978, 0.16238559048216278), (0.5133965847044396, 0.1724910290060434), (0.5123837853430446, 0.1824839666234244), (0.5113073982689569, 0.19260972531484918), (0.5101838879807972, 0.20267105739542704), (0.5089971915956811, 0.2127884725880221), (0.5077354074332255, 0.22303648323476166), (0.5064772997059829, 0.23274293946694702), (0.505135968407049, 0.2425720672408397), (0.5037744909162654, 0.25218755929500775), (0.502340183478718, 0.26187453857680526), (0.500789204001839, 0.272004056990892), (0.49919437872278605, 0.2820967459020566), (0.49755340494065897, 0.2921422920633047), (0.49582368341361727, 0.30237516359783323), (0.49421112885615637, 0.3115582643980614), (0.4925655065038384, 0.3206927992596466), (0.4908711492757497, 0.32979630941907917), (0.4891089752511264, 0.3390183222085125), (0.4872631520227969, 0.3483747130857823), (0.4853899983559906, 0.3576223959140362), (0.48346797983189377, 0.36682591159199557), (0.4814634844759115, 0.3762025640141872), (0.47914968818342946, 0.3868165904966582), (0.47674398195233264, 0.39746890356478165), (0.47428913546462825, 0.40810241221043575), (0.47175827476623466, 0.41876210891233506), (0.46917505098071083, 0.4292805449725414), (0.4664836788790888, 0.4400124878192158), (0.463777414519755, 0.45049725750310804), (0.4610057233571947, 0.4610057233571947)]

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
    text(" A1: {:.14f}".format(A1), 10, 56)
    text(" A2: {:.14f}".format(A2), 10, 72)
    text(" A3: {:.14f}".format(A3), 10, 88)
    text(" A4: {:.14f}".format(A4), 10, 104)
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
