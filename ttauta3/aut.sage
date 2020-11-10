var('a','b','c','d','x','y')
assume(a >= 0, b >= 0, c >= 0, d >= 0, x >= 0, y >=0)
automaton = {
    "abcx":[("sa", [a>=b ],{a:a+x-b, b:b, c:c, x:x}),
            ("sa", [b>=a       ],{a:x, b:a, c:c, x:x+b-a     }),
            ("sb", [b>=c, b-c>=x ],{a:a, b:a+b-c-x, c:c, x:x   }),
            ("sb", [b>=c, x>=b-c ],{a:a+x-b+c, b:a, c:c, x:b-c }),
            ("sc", [      ],{b:b,c:b+c+x, a:a, x:x     }),
            ("ta", [a>=x       ],{b:b, a:b+a-x, c:c, x:x     }),
            ("ta", [x>=a       ],{b:b+x-a, a:b, c:c, x:a     }),
            ("tb", [a>=b       ],{c:c, b:c+x, a:b, x:x+a-b   }),
            ("tb", [b>=a       ],{c:c, b:c+x+b-a, a:a, x:x   }),
            ("tc", [c>=b, c-b>=x ],{c:c-b-x, b:b, a:a, x:x     }),
            ("sb", [c>=b       ],{a:a+x, b:a, c:b, y:c-b     }),
            ("tc", [c>=b, x>=c-b ],{d:x-c+b, b:b, a:a, x:c-b   }),
            ("td", [      ],{a:a, d:a+x, b:b, x:x+c     }),
            ("tc", [b>=c       ],{d:x, b:c, a:a, y:b-c       }),
            ("sd", [      ],{c:c+x,d:c, b:b, y:a       })],
    "abcy" :[("tb", [a>=b       ],{c:c+y, b:c, a:b, x:a-b     }),
             ("sa", [a>=b, a-b>=y ],{a:a-b-y, b:b, c:c, y:y     }),
             ("sb", [b>=c       ],{a:a, b:a+y+b-c, c:c, y:y   }),
             ("sb", [c>=b       ],{a:a, b:a+y, c:b, y:y+c-b   }),
             ("sc", [c>=y       ],{b:b, c:b+c-y, a:a, y:y     }),
             ("sc", [y>=c       ],{b:b+y-c, c:b, a:a, y:c     }),
             ("ta", [      ],{b:b,a:b+a+y, c:c, y:y     }),
             ("tb", [b>=a, b-a>=y ],{c:c, b:c+b-a-y, a:a, y:y   }),
             ("tb", [b>=a, y>=b-a ],{c:c+y-b+a, b:c, a:a, y:b-a }),
             ("tc", [b>=c       ],{c:y, b:c, a:a, y:y+b-c     }),
             ("tc", [c>=b       ],{c:c+y-b, b:b, a:a, y:y     }),
             ("td", [      ],{a:a+y,d:a, b:b, x:c       }),
             ("sa", [b>=a       ],{d:y, b:a, c:c, x:b-a       }),
             ("sa", [a>=b, y>=a-b ],{d:y-a+b, b:b, c:c, y:a-b   }),
             ("sd", [      ],{c:c, d:c+y, b:b, y:y+a     })],
    "abdx" :[("sc", [      ],{b:b, c:b+x, a:a, x:x+d     }),
             ("sd", [d>=a, x>=d-a ],{c:x-d+a, a:a, b:b, x:d-a   }),
             ("sd", [a>=d       ],{c:x, a:d, b:b, y:a-d       }),
             ("sa", [a>=b       ],{d:d, a:d+x+a-b, b:b, x:x   }),
             ("sa", [b>=a       ],{d:d, a:d+x, b:a, x:x+b-a   }),
             ("sb", [b>=x       ],{a:a, b:a+b-x, d:d, x:x     }),
             ("sb", [x>=b       ],{a:a+x-b, b:a, d:d, x:b     }),
             ("sd", [d>=a, d-a>=x ],{d:d-a-x, a:a, b:b, x:x     }),
             ("ta", [a>=d, a-d>=x ],{b:b, a:b+a-d-x, d:d, x:x   }),
             ("ta", [a>=d, x>=a-d ],{b:b+x-a+d, a:b, d:d, x:a-d }),
             ("tb", [a>=b       ],{b:x, a:b, d:d, x:x+a-b     }),
             ("tb", [b>=a       ],{b:b+x-a, a:a, d:d, x:x     }),
             ("td", [      ],{a:a,d:a+d+x, b:b, x:x     }),
             ("ta", [d>=a       ],{b:b+x, a:b, d:a, y:d-a     }),
             ("tc", [      ],{d:d+x,c:d, a:a, y:b       })],
    "abdy" :[("sc", [      ],{b:b+y,c:b, a:a, x:d       }),
             ("sa", [b>=a       ],{d:d+y, a:d, b:a, x:b-a     }),
             ("sa", [a>=b, a-b>=y ],{d:d, a:d+a-b-y, b:b, y:y   }),
             ("sa", [a>=b, y>=a-b ],{d:d+y-a+b, a:d, b:b, y:a-b }),
             ("sb", [],{a:a,b:a+b+y, d:d, y:y     }),
             ("sd", [a>=d       ],{d:y, a:d, b:b, y:y+a-d     }),
             ("sd", [d>=a       ],{d:d+y-a, a:a, b:b, y:y     }),
             ("ta", [a>=d       ],{b:b, a:b+y+a-d, d:d, y:y   }),
             ("ta", [d>=a       ],{b:b, a:b+y, d:a, y:y+d-a   }),
             ("tb", [b>=a, b-a>=y ],{b:b-a-y, a:a, d:d, y:y     }),
             ("td", [d>=y       ],{a:a, d:a+d-y, b:b, y:y     }),
             ("td", [y>=d       ],{a:a+y-d, d:a, b:b, y:d     }),
             ("tb", [a>=b       ],{c:y, a:b, d:d, x:a-b       }),
             ("tb", [b>=a, y>=b-a ],{c:y-b+a, a:a, d:d, y:b-a   }),
             ("tc", [      ],{d:d, c:d+y, a:a, y:y+b     })],
    "acdx" :[("sb", [],{a:a+x,b:a, d:d, y:c       }),
             ("sa", [      ],{d:d,a:d+a+x, c:c, x:x     }),
             ("sc", [c>=d       ],{c:c+x-d, d:d, a:a, x:x     }),
             ("sc", [d>=c       ],{c:x, d:c, a:a, x:x+d-c     }),
             ("sd", [d>=a, d-a>=x ],{c:c, d:c+d-a-x, a:a, x:x   }),
             ("sd", [d>=a, x>=d-a ],{c:c+x-d+a, d:c, a:a, x:d-a }),
             ("ta", [a>=d, a-d>=x ],{a:a-d-x, d:d, c:c, x:x     }),
             ("tc", [c>=x       ],{d:d, c:d+c-x, a:a, x:x     }),
             ("tc", [x>=c       ],{d:d+x-c, c:d, a:a, x:c     }),
             ("td", [c>=d       ],{a:a, d:a+x, c:d, x:x+c-d   }),
             ("td", [d>=c       ],{a:a, d:a+x+d-c, c:c, x:x   }),
             ("sd", [a>=d       ],{c:c+x, d:c, a:d, y:a-d     }),
             ("ta", [a>=d, x>=a-d ],{b:x-a+d, d:d, c:c, x:a-d   }),
             ("tb", [      ],{c:c, b:c+x, d:d, x:x+a     }),
             ("ta", [d>=a       ],{b:x, d:a, c:c, y:d-a       })],
    "acdy" :[("sc", [d>=c       ],{b:y, d:c, a:a, x:d-c       }),
             ("sb", [],{a:a, b:a+y, d:d, y:y+c     }),
             ("sc", [c>=d, y>=c-d ],{b:y-c+d, d:d, a:a, y:c-d   }),
             ("td", [c>=d       ],{a:a+y, d:a, c:d, x:c-d     }),
             ("sa", [a>=y       ],{d:d, a:d+a-y, c:c, y:y     }),
             ("sa", [y>=a       ],{d:d+y-a, a:d, c:c, y:a     }),
             ("sc", [c>=d, c-d>=y ],{c:c-d-y, d:d, a:a, y:y     }),
             ("sd", [a>=d       ],{c:c, d:c+y, a:d, y:y+a-d   }),
             ("sd", [d>=a       ],{c:c, d:c+y+d-a, a:a, y:y   }),
             ("ta", [a>=d       ],{a:a+y-d, d:d, c:c, y:y     }),
             ("ta", [d>=a       ],{a:y, d:a, c:c, y:y+d-a     }),
             ("tc", [      ],{d:d,c:d+c+y, a:a, y:y     }),
             ("td", [d>=c, d-c>=y ],{a:a, d:a+d-c-y, c:c, y:y   }),
             ("td", [d>=c, y>=d-c ],{a:a+y-d+c, d:a, c:c, y:d-c }),
             ("tb", [      ],{c:c+y,b:c, d:d, x:a       })],
    "bcdx" :[("ta", [      ],{b:b+x,a:b, c:c, y:d       }),
             ("sa", [      ],{d:d, a:d+x, c:c, x:x+b     }),
             ("sb", [b>=c, x>=b-c ],{a:x-b+c, c:c, d:d, x:b-c   }),
             ("sb", [c>=b       ],{a:x, c:b, d:d, y:c-b       }),
             ("sb", [b>=c, b-c>=x ],{b:b-c-x, c:c, d:d, x:x     }),
             ("sc", [c>=d       ],{b:b, c:b+x+c-d, d:d, x:x   }),
             ("sc", [d>=c       ],{b:b, c:b+x, d:c, x:x+d-c   }),
             ("sd", [d>=x       ],{c:c, d:c+d-x, b:b, x:x     }),
             ("sd", [x>=d       ],{c:c+x-d, d:c, b:b, x:d     }),
             ("tb", [      ],{c:c,b:c+b+x, d:d, x:x     }),
             ("tc", [c>=b, c-b>=x ],{d:d, c:d+c-b-x, b:b, x:x   }),
             ("tc", [c>=b, x>=c-b ],{d:d+x-c+b, c:d, b:b, x:c-b }),
             ("td", [c>=d       ],{d:x, c:d, b:b, x:x+c-d     }),
             ("td", [d>=c       ],{d:d+x-c, c:c, b:b, x:x     }),
             ("tc", [b>=c       ],{d:d+x, c:d, b:c, y:b-c     })],
    "bcdy" :[("td", [c>=d       ],{a:y, c:d, b:b, x:c-d       }),
             ("ta", [      ],{b:b, a:b+y, c:c, y:y+d     }),
             ("td", [d>=c, y>=d-c ],{a:y-d+c, c:c, b:b, y:d-c   }),
             ("sa", [      ],{d:d+y,a:d, c:c, x:b       }),
             ("sc", [d>=c       ],{b:b+y, c:b, d:c, x:d-c     }),
             ("sb", [b>=c       ],{b:b+y-c, c:c, d:d, y:y     }),
             ("sb", [c>=b       ],{b:y, c:b, d:d, y:y+c-b     }),
             ("sc", [c>=d, c-d>=y ],{b:b, c:b+c-d-y, d:d, y:y   }),
             ("sc", [c>=d, y>=c-d ],{b:b+y-c+d, c:b, d:d, y:c-d }),
             ("sd", [      ],{c:c,d:c+d+y, b:b, y:y     }),
             ("tb", [b>=y       ],{c:c, b:c+b-y, d:d, y:y     }),
             ("tb", [y>=b       ],{c:c+y-b, b:c, d:d, y:b     }),
             ("tc", [b>=c       ],{d:d, c:d+y, b:c, y:y+b-c   }),
             ("tc", [c>=b       ],{d:d, c:d+y+c-b, b:b, y:y   }),
             ("td", [d>=c, d-c>=y ],{d:d-c-y, c:c, b:b, y:y     })],
}

# The following will try F(V)
# F is a conditional function F = (conditions, dictionary)
# V is data: dictionary
def tryApply(F, V):
    (conditions, function) = F
    if all([condition.substitute(V) for condition in conditions]):
        return {v: function[v].substitute(V) for v in function.keys()}
    else:
        return None

# [(target, mapname, dictionary)]
def nFoldCompose(n, source):
    if (n == 0):
        return [(source,"id", {var(x):var(x) for x in list(source)} )]
    previous = nFoldCompose(n-1, source)
    new = []
    for composition in previous:
        (target, mapname, dictionary) = composition
        newSource = target
        for edge in automaton[target]:
            (edgeName, edgeConditions, edgeDictionary) = edge
            newDictionary = tryApply((edgeConditions, edgeDictionary), dictionary)
            if newDictionary != None:
                newTarget = ''.join(sorted([str(x) for x in edgeDictionary.keys()]))
                new.append((newTarget, edgeName + "-o-" + mapname, newDictionary))
    return new

[('abcx', 'sc', {b: b, c: b + c + x, a: a, x: x}),
 ('abdx', 'td', {a: a, d: a + x, b: b, x: c + x}),
 ('bcdy', 'sd', {c: c + x, d: c, b: b, y: a}),

 ('abcy', 'sb-sc', {a: a + x, b: a, c: b, y: c + x}),
 ('abdy', 'ta-td', {b: b + c + x, a: b, d: a, y: x}),
 ('acdy', 'tc-td', {d: a + c + 2*x, c: a + x, a: a, y: b}),
 ('abcy', 'ta-sd', {b: b, a: a + b, c: c + x, y: a + c}),
 ('acdx', 'sa-sd', {d: a + c, a: c, c: c + x, x: b}),
 
 ('bcdy', 'sa-sb-sc', {d: c, b: a, c: b, y: x}),
 ('bcdy', 'tc-sd-sc',  {d: b + c + x, c: a + b + 2*c + 3*x, b: b, y: a}),
 ('abdy', 'tb-ta-td', {b: c, a: b, d: a, y: x}),
 ('acdy', 'sd-tc-td', {c: a + x, d: a + b + c + 3*x, a: a, y: b}),
 ('acdy', 'ta-tc-td', {a: b, d: a, c: a + x, y: b + c + 2*x}),
 ('abcx', 'tb-ta-sd', {c: a + 2*c + x, b: c + x, a: b, x: a})]
