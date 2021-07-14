from flask import Flask
from flask import request
from math import*
import cmath
app1 = Flask(__name__)
@app1.route("/")
def index0():
    return(        """
      <br><br> <center> <table>
     <tr>
    <td>  <a  href = "http://127.0.0.1:8080/hi" >Fourth degree solver</a> </td>
    
    <td>&nbsp;&nbsp;<td>
    <td>  <a  href = "http://127.0.0.1:8080/hi2" >Third degree solver</a> </td>
    <td>&nbsp;&nbsp;<td>
     <td>  <a  href = "http://127.0.0.1:8080/hi1" >Second degree solver</a> </td>
     <td>&nbsp;&nbsp;<td>
     <td>  <a  href = "http://127.0.0.1:8080/hi3" >Numerical solver</a> </td>

</table> <center>
        """)

@app1.route("/hi")
def index():
    a = request.args.get("a", "")
    b = request.args.get("b", "")
    c = request.args.get("c", "")
    d = request.args.get("d", "")
    e = request.args.get("e", "")
    if a:
        fahrenheit = fahrenheit_from(a, b, c , d ,e)
    else:
        fahrenheit = "", "", "", "","", "", "", "", ""
    return (
        """
        <form action="" method="get">
         If y = ax^4 + bx^3 + cx^2 + dx + e <br>
                 enter a: <input type="text" name="a"> <br>
                 enter b: <input type="text" name="b"> <br>
                 enter c: <input type="text" name="c"> <br>
                 enter d: <input type="text" name="d"> <br>
                 enter e: <input type="text" name="e">
                <input type="submit" value="Solve">
            </form>"""
        + "x = "
        + fahrenheit[0] + ", x = " + fahrenheit[1] + ", x = " + fahrenheit[2] + ', x = ' + fahrenheit[3]
        + "<br>Uitwerkingen: <br> " + fahrenheit[4] + fahrenheit[5] + fahrenheit[6] + fahrenheit[7] + fahrenheit[8])

def mypower(a,b):
    if a < 0:
        return -pow(abs(a), b)
    else:
        return pow(a,b)

def from_alpha_to_x(alpha,pin,qin, b):
    if round(alpha.imag,5) == 0:
        alpha = alpha.real

    gamma = 1 / 2 * (alpha**2 + pin + qin / alpha)
    beta = 1 / 2 * (alpha**2 + pin - qin / alpha)
    discriminant = round(alpha**2 - 4 * 1 * beta, 5)
    test = isinstance(discriminant, complex)
    if test == True:
        y1 = (-alpha + cmath.sqrt(discriminant)) / (2 * 1)
        y2 = (-alpha - cmath.sqrt(discriminant)) / (2 * 1)
        x1 = y1 - 1 / 4 * b
        x2 = y2 - 1 / 4 * b
        x1 = complex(round(x1.real,3), round(x1.imag, 3))
        x2 = complex(round(x2.real, 3), round(x2.imag, 3))
    else:
        if discriminant >= 0:
            y1 = (-alpha + sqrt(discriminant)) / (2 * 1)
            y2 = (-alpha - sqrt(discriminant)) / (2 * 1)
            x1 = round(y1 - 1 / 4 * b, 5)
            x2 = round(y2 - 1 / 4 * b, 5)
            x5 = None
        else:
            x1 = "No real answer"
            x2 = "No real answer"
            x5 = "No real answer"
    discriminant1 = alpha**2 - 4 * 1 * gamma
    if isinstance(discriminant1, complex) == True:
        y3 = (-alpha + cmath.sqrt(discriminant)) / (2 * 1)
        y4 = (-alpha - cmath.sqrt(discriminant)) / (2 * 1)
        x3 = y3 - 1 / 4 * b
        x4 = y4 - 1 / 4 * b
        x3 = complex(round(x3.real, 3), round(x3.imag, 3))
        x4 = complex(round(x4.real, 3), round(x4.imag, 3))
        x5 = complex(round(x1.real, 3), round(-x1.imag, 3))
    else:
        if discriminant1 >= 0:
            y3 = (alpha + sqrt(discriminant1)) / (2 * 1)
            y4 = (alpha - sqrt(discriminant1)) / (2 * 1)
            x3 = round(y3 - 1 / 4 * b, 5)
            x4 = round(y4 - 1 / 4 * b, 5)
            x5 = None
        else:
            x3 = "No real answer"
            x4 = "No real answer"
            x5 = "No real answer"
    return(x1, x2, x3, x4,x5, gamma, beta)

def cube_root(Uin):
    real_x = Uin.real
    imaginary_y = Uin.imag
    zpythagoras = sqrt(pow(imaginary_y, 2) + pow(real_x, 2))
    alpha = degrees(acos((real_x) / zpythagoras))
    beta = alpha / 3
    realpart = (cos(radians(beta))) * mypower(zpythagoras, 1 / 3)
    imaginarypart = (sin(radians(beta))) * mypower(zpythagoras, 1 / 3)
    cube_root = complex(realpart, imaginarypart)
    beta2 = 360 - 120 - 120 + beta
    realpart2 = (cos(radians(beta2))) * mypower(zpythagoras, 1 / 3)
    imaginarypart2 = (sin(radians(beta2))) * mypower(zpythagoras, 1 / 3)
    cube_root_2 = complex(realpart2, imaginarypart2)
    beta3 = 120 + 120 + beta
    realpart3 = (cos(radians(beta3))) * mypower(zpythagoras, 1 / 3)
    imaginarypart3 = (sin(radians(beta3))) * mypower(zpythagoras, 1 / 3)
    cube_root_3 = complex(realpart3, imaginarypart3)
    return (cube_root, cube_root_2, cube_root_3)

def print_values(xvalues):
    i = 0
    while i != len(xvalues):
        try:
            xvalues.remove("No real answer")
        except:
            pass
        print("x = " + str(xvalues[i]))
        i += 1

def a_lot_of_zeros(a3,a2,a1,a, pin, qin):
    Rt = 9 * a3 * a2 * a1 - 27 * pow(a3, 2) * a - 2 * pow(a2, 3)
    Rn = 54 * pow(a3, 3)
    R = Rt / Rn
    Qt = 3 * a3 * a1 - pow(a2, 2)
    Qn = 9 * pow(a3, 2)
    Q = Qt / Qn
    D = pow(Q, 3) + pow(R, 2)
    Sin = R + sqrt(abs(D))
    S = mypower(abs(Sin), 1 / 3)
    Tin = R - sqrt(abs(D))
    T = mypower(abs(Tin), 1 / 3)
    realans = sqrt((S + T) + (- a2 / (3 * a3)))
    if realans == 0:
        beta = 1 / 2 * (pow(realans, 2) + pin - qin)
        gamma = 1 / 2 * (pow(realans, 2) + pin + qin)
    else:
        beta = 1 / 2 * (pow(realans, 2) + pin - qin / realans)
        gamma = 1 / 2 * (pow(realans, 2) + pin + qin / realans)
    discriminant1 = pow(realans, 2) - 4 * 1 * gamma
    y3 = (realans + sqrt(discriminant1)) / (2 * 1)
    y4 = (realans - sqrt(discriminant1)) / (2 * 1)
    x3 = round(y3 - 1 / 4 * b, 5)
    x4 = round(y4 - 1 / 4 * b, 5)
    discriminant = pow(realans, 2) - 4 * 1 * beta
    y1 = (-realans + sqrt(discriminant)) / (2 * 1)
    y2 = (-realans - sqrt(discriminant)) / (2 * 1)
    x1 = round(y1 - 1 / 4 * b, 5)
    x2 = round(y2 - 1 / 4 * b, 5)
    return(x3, x4, x1, x2)

def calculations(a0, b, c, d, e, pin, qin, rin, a3, a2, a1, a, p, q,
                 Udoublein, Uin, Uin0, cuberoots,
                 ans, ans1, ans2, alpha,
                 alpha1, alpha2, gamma0,
                 beta0, gamma1, beta1,
                 gamma2, beta2, check, xvalues, xvalues1, xvalues2):
    uitwerkingpart1 = ("The equation is {a}x^4 + {b}x^3 + {c}x^2 + {d}x + {e} <br>"
                       "Divide all coefficients by a to reduce the equation to x^4 + b'x^3 + c'x^2 + d'x + e'."
                       "<br>b/a = {ba}, c/a = {ca}, d/a = {da}, e/a = {ea}<br>"
                       "so: x^4 + {ba}x^3 + {ca}x^2 + {da}^x + {ea}<br>"
                       "Using the transformation x = y - 1/4b, the equation can be rewritten without the third power in the form y^4 + py^2 + qy +r."
                       "<br> with p = c' - 3/8 * b'^2 = {ca} - 3/8 * {ba}^2 = {p},"
                       "<br>      q = 1/8 * b'^3 - 1/2 * b'*c' + d = 1/8 * {ba}^3 - 1/2 * {ba} * {ca} + {da} = {q}"
                       "<br> and  r = - 3/256 * b'^4 + 1/16 * b'^2 - 1/4b'*d' + e = -3/256 * {ba}^4 + 1/16 * {ba}^2 - 1/4 * {ba} * {da} + {ea} = {r}  "
                       "<br>This equation can be rewritten as a product of two quadratic equations"
                       "<br>y^4 + {p}y^2 + {q}y + {r} = (y^2 + αy + β)(y^2 - αy + γ) = y^4 + (γ + β - α^2)y^2 + α(γ - β)y + β*y"
                       "<br>From this you can conclude that p = γ + β - α^2, q = α(γ - β) r = β*y"
                       "<br>So: p + α^2 = γ + β, q/α = γ - β which you can rewrite as:"
                       "<br>4r = (p + α^2)^2 - (q/α)^2"
                       "<br>α^2 = A"
                       "<br>A^3 + (2p)A^2 + (p^2-4r)A - q^2"
                       "<br>= 0"
                       "<br>This is a third degree equation in α^2 (A)"
                       "<br>So in this case the equation is as follows:"
                       "<br>A^3 + {a2}A^2 + {a1}A - {a0} = 0"
                       "<br>This equation is solvable using the Cardano formula which in this case is:"
                       "<br>x = u - p/3u where u = 3√(-q/2 + √((q^2)/4 + (p^3)/27))"
                       "<br>First, using the transformation x = y + a/3, the equation can be rewritten without the second power in the form y^3 + py + q."
                       "<br>p = c - (b^2)/3 = {a1} - {a2}^2/3 = {p3}"
                       "<br>q = d + (2* b^3)/27 - (c*b)/3 = {a0} + (2* {a2}^3)/27 = {q3} "
                       "<br>In this cube root, the square root can be calculated first: 3√(-q/2 + √((q^2)/4 + (p^3)/27))"
                       "<br> √((q^2)/4 + (p^3)/27)) = √(({q3}^2)/4 + ({p3}^3)/27)) = √{Udoublein}".format(a=a0, b=b, c=c,
                                                                                                        d=d, e=e,
                                                                                                        ba=b / a0,
                                                                                                        ca=c / a0,
                                                                                                        da=d / a0,
                                                                                                        ea=e / a0,
                                                                                                        p=round(pin, 3),
                                                                                                        q=round(qin, 3),
                                                                                                        r=round(rin, 3),
                                                                                                        a2=round(a2, 3),
                                                                                                        a1=round(a1, 3),
                                                                                                        a0=round(a, 3),
                                                                                                        p3=round(p, 3),
                                                                                                        q3=round(q, 3),
                                                                                                        Udoublein=round(
                                                                                                            Udoublein,
                                                                                                            3)))
    if check == 1:
        uitwerkingpart2 = ("<br>In this case, the square root is negative so it's an complex number:"
                           "<br> Then the -q/2 is added to it:"
                           "<br> √{Udoublein} - {q3}/2 = {Uin} "
                           "<br> Then you take the cube root of this complex number which gives ="
                           "<br> y = 3√({Uin}) - {p3}/3 x 3√({Uin}) = {cuberoot0} - {p3}/(3*{cuberoot0} = {ans}"
                           "<br> v y = 3√({Uin}) - {p3}/3 x 3√({Uin}) = {cuberoot1} - {p3}/(3*{cuberoot1} = {ans1} "
                           "<br> v y = y = 3√({Uin}) - {p3}/3 x 3√({Uin}) = {cuberoot2} - {p3}/(3*{cuberoot2} = {ans2} "
                           "<br> A is y + a/3"
                           "<br> Dat geeft:"
                           "<br> α = √({ans} - {a2}/3) = {alpha} "
                           "<br> v α = √{ans1} - {a2}/3) = {alpha1} "
                           "<br> v α = √({ans2} - {a2}/3) =  {alpha2}"
                           "".format(a2=round(a2, 3), p3=round(p, 3), q3=round(q, 3), Udoublein=round(Udoublein, 3),
                                     Uin=Uin0,
                                     cuberoot0=complex(round(cuberoots[0].real, 3), round(cuberoots[0].imag, 3)),
                                     cuberoot1=complex(round(cuberoots[1].real, 3), round(cuberoots[1].imag, 3)),
                                     cuberoot2=complex(round(cuberoots[2].real, 3), round(cuberoots[2].imag, 3))
                                     , ans=round(ans.real, 3), ans1=round(ans1.real, 3), ans2=round(ans2.real, 3),
                                     alpha=round(alpha.real, 3), alpha1=round(alpha1.real, 3),
                                     alpha2=round(alpha2.real, 3)))
        if round(alpha.real,3) != 0:
            uitwerkingpart3 = ("From p + α^2 = γ + β, q/α = γ - β, you can conclude the following:"
                  "<br> For this alpha value:"
                  "<br>γ = 1/2(α^2 + p + q/α) = 1/2({alpha}^2 + {p} + {q}/{alpha} = {gamma}"
                  "<br>β = 1/2(α^2 + p - q/α) = 1/2({alpha}^2 + {p} - {q}/{alpha} = {beta}"
                  "<br> y^2 + αy + β = 0 and y^2 - αy + γ are simple quadratic equations."
                  "<br> In this case it is: y^2 + {alpha}y + {beta} (remember x = y - 1/4b)"
                  "<br> x = (-{alpha} + √({alpha}^2 - 4 * {beta}))/2 - 1/4 * {ba} = {x1} "
                  "<br> x = (-{alpha} - √({alpha}^2 - 4 * {beta}))/2 - 1/4 * {ba} = {x2}"
                  "<br> and y^2 - {alpha}y + {gamma} "
                  "<br> x = ({alpha} + √({alpha}^2 - 4 * {gamma}))/2 - 1/4 * {ba} = {x3}"
                  "<br> x = ({alpha} - √({alpha}^2 - 4 * {gamma}))/2 - 1/4 * {ba} = {x4}"
                  " ".format(x3 = xvalues[2], x4 = xvalues[3] ,x1 = xvalues[0], x2 = xvalues[1],ba= b/a0, alpha = round(alpha.real,3), p = pin, q = qin, beta = beta0, gamma = gamma0))
        else:
            uitwerkingpart3 = ""
        if round(alpha1.real,3) != 0:
            uitwerkingpart4 = ("From p + α^2 = γ + β, q/α = γ - β, you can conclude the following:"
                  "<br> For this alpha value:"
                  "<br>γ = 1/2(α^2 + p + q/α) = 1/2({alpha}^2 + {p} + {q}/{alpha} = {gamma}"
                  "<br>β = 1/2(α^2 + p - q/α) = 1/2({alpha}^2 + {p} - {q}/{alpha} = {beta}"
                  "<br> y^2 + αy + β = 0 and y^2 - αy + γ are simple quadratic equations."
                  "<br> In this case it is: y^2 + {alpha}y + {beta} (remember x = y - 1/4b)"
                  "<br> x = (-{alpha} + √({alpha}^2 - 4 * {beta}))/2 - 1/4 * {ba} = {x1} "
                  "<br> x = (-{alpha} - √({alpha}^2 - 4 * {beta}))/2 - 1/4 * {ba} = {x2}"
                  "<br> and y^2 - {alpha}y + {gamma} "
                  "<br> x = ({alpha} + √({alpha}^2 - 4 * {gamma}))/2 - 1/4 * {ba} = {x3}"
                  "<br> x = ({alpha} - √({alpha}^2 - 4 * {gamma}))/2 - 1/4 * {ba} = {x4}"
                  " ".format(x3=xvalues1[2], x4=xvalues1[3], x1=xvalues1[0], x2=xvalues1[1], ba=b / a0,
                             alpha=round(alpha1.real, 3), p=pin, q=qin, beta=beta1, gamma=gamma1))
        else:
            uitwerkingpart4 = ""
        if round(alpha2.real,3) != 0:
            uitwerkingpart5 = (
                  "<br><br><br> For this alpha value:"
                  "<br>γ = 1/2(α^2 + p + q/α) = 1/2({alpha}^2 + {p} + {q}/{alpha} = {gamma}"
                  "<br>β = 1/2(α^2 + p - q/α) = 1/2({alpha}^2 + {p} - {q}/{alpha} = {beta}"
                  "<br> y^2 + αy + β = 0 and y^2 - αy + γ are simple quadratic equations."
                  "<br> In this case it is: y^2 + {alpha}y + {beta} (remember x = y - 1/4b)"
                  "<br> x = (-{alpha} + √({alpha}^2 - 4 * {beta}))/2 - 1/4 * {ba} = {x1} "
                  "<br> x = (-{alpha} - √({alpha}^2 - 4 * {beta}))/2 - 1/4 * {ba} = {x2}"
                  "<br> and y^2 - {alpha}y + {gamma} "
                  "<br> x = ({alpha} + √({alpha}^2 - 4 * {gamma}))/2 - 1/4 * {ba} = {x3}"
                  "<br> x = ({alpha} - √({alpha}^2 - 4 * {gamma}))/2 - 1/4 * {ba} = {x4}"
                  " ".format(x3=xvalues2[2], x4=xvalues2[3], x1=xvalues2[0], x2=xvalues2[1], ba=b / a0,
                             alpha=round(alpha2.real, 3), p=pin, q=qin, beta=beta2, gamma=gamma2))
        else:
            uitwerkingpart5 = ""

    if check == 3:
        uitwerkingpart2 = ("<br>In this case, the square root is 0 so this is a special case:"
              "<br> first the -q/2 is added to it:"
              "<br> √{Udoublein} - {q3}/2 = {Uin}"
              "<br> In this case:"
              "<br> (regular) y = 3√{Uin} - {a2}/(3√{Uin} * 3)  = {ans}"
              "<br> (special) y = -3√{Uin} = {ans1}"
              "<br>A is y + a/3"
              "<br>Dat geeft:"
              "<br> α = √({ans} - {a2}/3) = {alpha} "
              "<br> v α = √{ans1} - {a2}/3) = {alpha1} "
              "".format(a2= round(a2,3),Udoublein = round(Udoublein,3), q3 = round(q,3), Uin = round(Uin,3)
                        , ans = round(ans,3), ans1 = round(ans1,3), alpha = round(alpha,3), alpha1= round(alpha1,3) ))
        if round(alpha, 3) != 0:
            uitwerkingpart3 = ("From p + α^2 = γ + β, q/α = γ - β, you can conclude the following:"
                  "<br> For this alpha value:"
                  "<br>γ = 1/2(α^2 + p + q/α) = 1/2({alpha}^2 + {p} + {q}/{alpha} = {gamma}"
                  "<br>β = 1/2(α^2 + p - q/α) = 1/2({alpha}^2 + {p} - {q}/{alpha} = {beta}"
                  "<br> y^2 + αy + β = 0 and y^2 - αy + γ are simple quadratic equations."
                  "<br> In this case it is: y^2 + {alpha}y + {beta} (remember x = y - 1/4b)"
                  "<br> x = (-{alpha} + √({alpha}^2 - 4 * {beta}))/2 - 1/4 * {ba} = {x1} "
                  "<br> x = (-{alpha} - √({alpha}^2 - 4 * {beta}))/2 - 1/4 * {ba} = {x2}"
                  "<br> and y^2 - {alpha}y + {gamma} "
                  "<br> x = ({alpha} + √({alpha}^2 - 4 * {gamma}))/2 - 1/4 * {ba} = {x3}"
                  "<br> x = ({alpha} - √({alpha}^2 - 4 * {gamma}))/2 - 1/4 * {ba} = {x4}"
                  " ".format(x3=xvalues[2], x4=xvalues[3], x1=xvalues[0], x2=xvalues[1], ba=b / a0,
                             alpha=round(alpha.real, 3), p=pin, q=qin, beta=beta0, gamma=gamma0))
        else: uitwerkingpart3 = ""
        if round(alpha1.real,3) != 0:
            uitwerkingpart4 = ("From p + α^2 = γ + β, q/α = γ - β, you can conclude the following:"
                  "<br>For this alpha value:"
                  "<br>γ = 1/2(α^2 + p + q/α) = 1/2({alpha}^2 + {p} + {q}/{alpha} = {gamma}"
                  "<br>β = 1/2(α^2 + p - q/α) = 1/2({alpha}^2 + {p} - {q}/{alpha} = {beta}"
                  "<br>y^2 + αy + β = 0 and y^2 - αy + γ are simple quadratic equations."
                  "<br>In this case it is: y^2 + {alpha}y + {beta} (remember x = y - 1/4b)"
                  "<br>x = (-{alpha} + √({alpha}^2 - 4 * {beta}))/2 - 1/4 * {ba} = {x1} "
                  "<br>x = (-{alpha} - √({alpha}^2 - 4 * {beta}))/2 - 1/4 * {ba} = {x2}"
                  "<br>and y^2 - {alpha}y + {gamma} "
                  "<br>x = ({alpha} + √({alpha}^2 - 4 * {gamma}))/2 - 1/4 * {ba} = {x3}"
                  "<br>x = ({alpha} - √({alpha}^2 - 4 * {gamma}))/2 - 1/4 * {ba} = {x4}"
                  " ".format(x3=xvalues1[2], x4=xvalues1[3], x1=xvalues1[0], x2=xvalues1[1], ba=b / a0,
                             alpha=round(alpha1.real, 3), p=pin, q=qin, beta=beta1, gamma=gamma1))
        else: uitwerkingpart4 = ""
    uitwerkingpart5 = ""
    if check == 2:
        uitwerkingpart2 = ("<br>First the -q/2 is added to it:"
              "<br>√{Udoublein} - {q3}/2 = {Uin}"
              "<br>In this case:"
              "<br>y = 3√{Uin} - {a2}/(3√{Uin} * 3)  = {ans}"
              "<br>A is y + a/3"
              "<br>Dat geeft:"
              "<br>α = √({ans} - {a2}/3) = {alpha} ".format(Udoublein = Udoublein, q3 = q, Uin = Uin, ans = ans, a2 = a2, alpha = alpha ))
        uitwerkingpart3 = ("From p + α^2 = γ + β, q/α = γ - β, you can conclude the following:"
              "<br> For this alpha value:"
              "<br>γ = 1/2(α^2 + p + q/α) = 1/2({alpha}^2 + {p} + {q}/{alpha} = {gamma}"
              "<br>β = 1/2(α^2 + p - q/α) = 1/2({alpha}^2 + {p} - {q}/{alpha} = {beta}"
              "<br> y^2 + αy + β = 0 and y^2 - αy + γ are simple quadratic equations."
              "<br> In this case it is: y^2 + {alpha}y + {beta} (remember x = y - 1/4b)"
              "<br> x = (-{alpha} + √({alpha}^2 - 4 * {beta}))/2 - 1/4 * {ba} = {x1} "
              "<br> x = (-{alpha} - √({alpha}^2 - 4 * {beta}))/2 - 1/4 * {ba} = {x2}"
              "<br> and y^2 - {alpha}y + {gamma} "
              "<br> x = ({alpha} + √({alpha}^2 - 4 * {gamma}))/2 - 1/4 * {ba} = {x3}"
              "<br> x = ({alpha} - √({alpha}^2 - 4 * {gamma}))/2 - 1/4 * {ba} = {x4}".format(x3 = xvalues[2], x4 = xvalues[3] ,x1 = xvalues[0], x2 = xvalues[1],ba= b/a0, alpha = round(alpha.real,3), p = pin, q = qin, beta = beta0, gamma = gamma0))
        uitwerkingpart4 = ""
        uitwerkingpart5 = ""
    return uitwerkingpart1, uitwerkingpart2, uitwerkingpart3, uitwerkingpart4, uitwerkingpart5
def fahrenheit_from(a,b,c,d,e):
    a0 = float(a)
    b = (float(b) / float(a))
    c = (float(c) / float(a))
    d = (float(d) / float(a))
    e = (float(e) / float(a))
    pin = c - (3/8) * pow(b, 2)
    qin = (1/8) * pow(b, 3) - (1/2) * b * c + d
    rin = (-3/256) * pow(b, 4) + (1/16) * pow(b, 2) * c - (1/4) * b * d + e

    a3 = 1
    a2 = 2 * pin
    a1 = pow(pin, 2) -4 * rin
    a = -pow(qin, 2)
    p = a1 - (pow(a2, 2))/(3)
    q = a + (2 * pow(a2, 3))/27 - (a1 * a2/3)



    Udoublein = round(pow(q,2)/4 + pow(p, 3 )/27,5)

    if round(Udoublein) < 0:
        check = 1
        Uin = -q / 2 + cmath.sqrt(Udoublein)
        Uin0 = complex(round(Uin.real,3), round(Uin.imag,4))
        cuberoots = cube_root(Uin)


        ans = cuberoots[0] - p / (cuberoots[0] * 3)
        alpha = cmath.sqrt(round((ans - a2/3).real,5))
        if alpha == 0 or round(alpha.imag,3) != 0 :
            xvalues = []
            gamma0 = None
            beta0 = None
        else:
            xvalues = list(from_alpha_to_x(alpha, pin, qin, b))
            beta0 = xvalues[6]
            gamma0 = xvalues[5]
            del xvalues[5:7]
            v = 1


        ans1 = cuberoots[1] - p / (cuberoots[1] * 3)
        alpha1 = cmath.sqrt(round((ans1 - a2 / 3).real, 5))
        if alpha1 == 0 or round(alpha1.imag,3) != 0 :
            xvalues1 = []
            gamma1 = []
            beta1 = []
        else:
            xvalues1 = list(from_alpha_to_x(alpha1, pin, qin, b))
            gamma1 = xvalues1[5]
            beta1 = xvalues1[6]
            del xvalues1[5:7]
            v += 1
        ans2 = cuberoots[2] - p / (cuberoots[2] * 3)
        alpha2 = cmath.sqrt(round((ans2 - a2/3).real,5))
        if alpha2 == 0 or round(alpha2.imag,3) != 0 :
            xvalues2 = []
            gamma2 = []
            beta2 = []
        else:
            xvalues2 = list(from_alpha_to_x(alpha2, pin, qin, b))
            gamma2 = xvalues2[5]
            beta2 = xvalues2[6]
            del xvalues2[5:7]
            v +=1
        xvalues = list(set(xvalues + xvalues1 + xvalues2))
        uitwerkingen = calculations(a0, b, c, d, e, pin, qin, rin, a3, a2, a1, a, p, q,
                                    Udoublein, Uin, Uin0, cuberoots,
                                    ans, ans1, ans2, alpha,
                                    alpha1, alpha2, gamma0,
                                    beta0, gamma1, beta1,
                                    gamma2, beta2, check, xvalues, xvalues1, xvalues2)
        try:
            xvalues.remove("No real answer")
        except:
            pass
        while len(xvalues) < 5:
            xvalues.append("")

        return (str(xvalues[0]), str(xvalues[1]), str(xvalues[2]), str(xvalues[3]),
                str(uitwerkingen[0]), str(uitwerkingen[1]), str(uitwerkingen[2]), str(uitwerkingen[3]), str(uitwerkingen[4]))
    elif round(Udoublein, 5) > 0:
        check = 2
        Uin = -q / 2 + sqrt(Udoublein)
        U = mypower(Uin, 1 / 3)
        ans3 = U - p / (U * 3)
        realans = ans3 - a2 / 3
        alpha = round(sqrt(realans),3)
        xvalues = list((from_alpha_to_x(alpha, pin, qin, b)))
        gamma0 = xvalues[5]
        beta0 = xvalues[6]
        del xvalues[5:7]


        uitwerkingen = calculations(a0, b, c, d, e, pin, qin, rin, a3, a2, a1, a, p, q,
                 Udoublein, Uin, None, None,
                 ans3, None, None, alpha,
                 None, None, gamma0,
                 beta0, None, None,
                 None, None, check, xvalues, None, None)
        xvaluesfinal = list(set(xvalues))
        try:
            xvaluesfinal.remove("No real answer")
        except:
            pass
        while len(xvaluesfinal) < 5:
            xvaluesfinal.append("")

        return (str(xvaluesfinal[0]), str(xvaluesfinal[1]), str(xvaluesfinal[2]), str(xvaluesfinal[3]),
                str(uitwerkingen[0]), str(uitwerkingen[1]), str(uitwerkingen[2]), str(uitwerkingen[3]),
                str(uitwerkingen[4]))
    elif round(Udoublein, 5) == 0 and p != 0 and q != 0:
        check = 3
        Uin = -q / 2 + sqrt(Udoublein)
        U = mypower(Uin, 1 / 3)
        ans3 = U - p / (U * 3)
        ans4 = -U
        realans = round(ans3 - a2 / 3, 5)
        realans = sqrt(realans)
        if realans == 0:
            xvalues = []
            gamma0 = None
            beta0 = None
        else:
            xvalues = list(from_alpha_to_x(realans, pin, qin, b))
            gamma0 = xvalues[5]
            beta0 = xvalues[6]
            del xvalues[5:7]
        realans1 = round(ans4 - a2 / 3, 5)
        realans1 = sqrt(realans1)
        if realans1 == 0:
            xvalues1 = []
            gamma1 = None
            beta1 = None
        else:
            xvalues1 = list(from_alpha_to_x(realans1, pin, qin, b))
            beta1 = xvalues1[6]
            gamma1 = xvalues1[5]
            del xvalues1[5:7]
        xvaluesfinal = list(set(xvalues + xvalues1))
        try:
            xvaluesfinal.remove("No real answer")
        except:
            pass
        while len(xvaluesfinal) < 5:
            xvaluesfinal.append("")
        uitwerkingen = calculations(a0, b, c, d, e, pin,qin,rin,a3,a2,a1,a,p,q,
                     Udoublein, Uin, None, None,
                     ans3, ans4, None, realans,
                     realans1, None, gamma0,
                     beta0, gamma1, beta1
                     ,None,None, check, xvalues, xvalues1, None)
        return(str(xvaluesfinal[0]), str(xvaluesfinal[1]), str(xvaluesfinal[2]), str(xvaluesfinal[3]),
               str(uitwerkingen[0]), str(uitwerkingen[1]), str(uitwerkingen[2]), str(uitwerkingen[3]), str(uitwerkingen[4]))

    elif p == q == 0:
        check = 4
        xvalues = list(set(a_lot_of_zeros(a3,a2,a1,a,pin,qin)))
        while len(xvalues) < 5:
            xvalues.append("")
        uitwerkingen = "geenuitwerkingen"
        return (xvalues[0], xvalues[1], xvalues[2], xvalues[3], uitwerkingen, '','','', '')


@app1.route("/hi1")
def index1():
    aa = request.args.get("a", "")
    bb = request.args.get("b", "")
    cc = request.args.get("c", "")
    if aa:
        solution = solver(aa, bb, cc)
    else:
        solution = "", ""
    return (
        """<form action="" method="get">
         If y = ax^2 + bx + c <br>
                 enter a: <input type="text" name="a"> <br>
                 enter b: <input type="text" name="b"> <br>
                 enter c: <input type="text" name="c">
                <input type="submit" value="Solve">
            </form>"""
        + "x = "
        + solution[0] + " and x = " + solution[1]
    )

def solver(m,n,o):
    a = float(m)
    b = float(n)
    c = float(o)
    """Convert Celsius to Fahrenheit degrees."""
    discriminant = pow(b, 2) - 4 * a * c
    if discriminant < 0:
        x = (-b + cmath.sqrt(discriminant)) / (2 * a)
        x = complex(round(x.real, 3), round(x.imag,3))
        x1 = (-b - cmath.sqrt(discriminant)) / (2 * a)
        x1 = complex(round(x1.real, 3), round(x1.imag, 3))
        return str(x), str(x1)
    elif discriminant > 0:
        x = round((-b + sqrt(discriminant)) / (2 * a), 3)
        x1 = round((-b - sqrt(discriminant)) / (2 * a), 3)
        return str(x), str(x1)
    elif discriminant == 0:

        x = round((-b + sqrt(discriminant)) / (2 * a),3)
        return str(x), str(x)


@app1.route("/hi2")
def index2():
    aaa = request.args.get("a", "")
    bbb = request.args.get("b", "")
    ccc = request.args.get("c", "")
    ddd = request.args.get("d", "")
    if aaa:
        oplossing = derdegraadsoplosser(aaa, bbb, ccc, ddd)
    else:
        oplossing = "", "", ""
    return (
            """<form action="" method="get">
             If y = ax^3 + bx^2 + cx + d <br>
                     enter a: <input type="text" name="a"> <br>
                     enter b: <input type="text" name="b"> <br>
                     enter c: <input type="text" name="c"> <br>
                     enter d: <input type="text" name="d">
                    <input type="submit" value="Solve">
                </form>"""
            + "x = "
            + oplossing[0] + ", x = " + oplossing[1] + " and x = " + oplossing[1]
    )
def derdegraadsoplosser(a,b,c,d):
    a = float(a)
    b = (float(b) / float(a))
    c = (float(c) / float(a))
    d = (float(d) / float(a))
    p = c - (pow(b, 2)) / (3)
    q = d + (2 * pow(b, 3)) / 27 - (c * b / 3)

    Udoublein = pow(q, 2) / (4) + pow(p, 3) / (27)
    if round(Udoublein) < 0:
        Uin = -q / 2 + cmath.sqrt(Udoublein)

        real_x = Uin.real
        imaginary_y = Uin.imag
        zpythagoras = sqrt(pow(imaginary_y, 2) + pow(real_x, 2))
        alpha = degrees(acos((real_x) / zpythagoras))
        beta = alpha / 3
        realpart = (cos(radians(beta))) * mypower(zpythagoras, 1 / 3)
        imaginarypart = (sin(radians(beta))) * mypower(zpythagoras, 1 / 3)
        cube_root = complex(realpart, imaginarypart)
        beta2 = 360 - 120 - 120 + beta
        realpart2 = (cos(radians(beta2))) * mypower(zpythagoras, 1 / 3)
        imaginarypart2 = (sin(radians(beta2))) * mypower(zpythagoras, 1 / 3)
        cube_root_2 = complex(realpart2, imaginarypart2)
        beta3 = 120 + 120 + beta
        realpart3 = (cos(radians(beta3))) * mypower(zpythagoras, 1 / 3)
        imaginarypart3 = (sin(radians(beta3))) * mypower(zpythagoras, 1 / 3)
        cube_root_3 = complex(realpart3, imaginarypart3)

        ans = cube_root - p / (cube_root * 3)
        realans = ans - b / 3
        realans = realans.real

        ans1 = (cube_root_2 - p / (cube_root_2 * 3))
        realans1 = ans1 - b / 3
        realans1 = realans1.real

        ans2 = cube_root_3 - p / (cube_root_3 * 3)
        realans2 = ans2 - b / 3
        realans2 = realan2.real
        return(str(realans), str(realans1), str(realans2))
    elif round(Udoublein, 5) > 0:
        Uin = -q / 2 + sqrt(Udoublein)
        U = mypower(Uin, 1 / 3)
        ans3 = U - p / (U * 3)
        realans3 = ans3 - b / 3
        return(str(realans3), "", "" )

    elif round(Udoublein, 5) == 0 and p != 0 and q != 0:
        print(p, q)
        Uin = -q / 2 + sqrt(Udoublein)
        U = mypower(Uin, 1 / 3)
        ans3 = U - p / (U * 3)
        ans4 = -U
        realans3 = ans3 - b / 3
        realans4 = ans4 - b / 3
        return(str(realans3), str(realans4), "")

    elif p == q == 0:
        print("aaaa")
        Rt = 9 * a * b * c - 27 * pow(a, 2) * d - 2 * pow(b, 3)
        Rn = 54 * pow(a, 3)
        R = Rt / Rn
        Qt = 3 * a * c - pow(b, 2)
        Qn = 9 * pow(a, 2)
        Q = Qt / Qn
        D = pow(Q, 3) + pow(R, 2)
        Sin = R + sqrt(abs(D))
        S = mypower(abs(Sin), 1 / 3)
        Tin = R - sqrt(abs(D))
        T = mypower(abs(Tin), 1 / 3)

        realans5 = (S + T) + (- b / (3 * a))
        return(str(realans5), "", "")

@app1.route("/hi3")
def index3():
    return ("Not available yet  ")

if __name__ == "__main__":
    app1.run(host="127.0.0.1", port=8080, debug=True)
