import turtle as t

def recursionPrint(text, level):
    print(level * "\t\t" + text)

def treeSimpleVerbose(length, angle, minlength=10, level=0):
    recursionPrint("level: {} | length: {} | angle: {}".format(level, length, angle), level)

    if length < minlength:
        recursionPrint("length < {}, returning\n".format(minlength), level)
        return

    recursionPrint("forward({})".format(length), level)
    t.forward(length)
    recursionPrint("left({})".format(angle), level)
    t.left(angle)

    recursionPrint("calling first tree({}, {}, {})\n".format(length * 0.75, angle, level+1), level)
    treeSimpleVerbose(length * 0.75, angle, minlength=minlength,  level=level+1)

    recursionPrint("right({})".format(2 * angle), level)
    t.right(2 * angle)

    recursionPrint("calling second tree({}, {}, {})\n".format(length * 0.75, angle, level+1), level)
    treeSimpleVerbose(length * 0.75, angle, minlength=minlength, level=level+1)

    recursionPrint("left({})".format(angle), level)
    t.left(angle)
    recursionPrint("back({})".format(length), level)
    t.back(length)

def tree(length, angle, color, width):
    t.color(color)
    t.width(width)
    t.pendown()

    if length < 10:
        return

    t.forward(length)
    t.left(angle)

    tree(length*0.75, angle, (color[0], color[1] + 3, color[2] + 1), width * 1.03)

    t.right(2 * angle)

    tree(length*0.75, angle, (color[0], color[1] + 3, color[2] + 1), width * 1.03)

    t.left(angle)
    t.color(gcolor)
    t.penup()
    t.back(length)

glength = 200
gangle = 20
minlength = 50
gcolor = (168, 106, 39)
width = 1.5

t.speed(speed="normal")

t.left(90)
t.penup()
t.back(400)
t.pendown()

treeSimpleVerbose(glength, gangle, minlength=minlength) # verbose drawing
#tree(glength, gangle, gcolor, width)

t.mainloop()