import turtle

smallTree = ("ancestor",
                ("ancestor",
                    ("D", (), ()),
                    ("E", (), ())
                ),
                ("ancestor",
                    ("F", (), ()),
                    ("G", (), ())
                )
            )

smallTree2 = ("ancestor",
                ("ancestor",
                    ("ancestor",
                     ("ancestor",
                        ("Z", (), ()),
                        ("X", (), ())
                    ),
                    ("Y", (), ()))
                 ,
                    ("ancestor",
                     ("E", (), ()),
                     ("D", (), ()))
                ),
                ("ancestor",
                    ("F", (), ()),
                    ("ancestor",
                     ("A", (), ()),
                     ("B", (), ())
                    )
                )
            )

def drawTree(tree, angle, length, width):
    turtle.width(width)

    if tree[0] == "ancestor":
        # left branch
        turtle.left(angle)
        turtle.forward(length)
        turtle.right(angle)
        drawTree(tree[1], angle - 0.2 * angle, length - 0.2 * length, width - 0.3 * width)
        turtle.width(width)
        turtle.left(angle)
        turtle.backward(length)
        turtle.right(angle)
        
        # right branch
        turtle.right(angle)
        turtle.forward(length)
        turtle.left(angle)
        drawTree(tree[2], angle - 0.2 * angle, length - 0.2 * length, width - 0.3 * width)
        turtle.width(width)
        turtle.right(angle)
        turtle.backward(length)
        turtle.left(angle)
    else:
        # draw the ending node
        turtle.pencolor("red")
        turtle.write(tree[0], font=("Monospace", 14, "bold"))
        turtle.pencolor("black")

turtle.speed(speed="normal")
drawTree(smallTree2, 45, 100, 5)
turtle.mainloop()