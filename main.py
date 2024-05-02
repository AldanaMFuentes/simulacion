import math
import random
import matplotlib.pyplot as plt
import scipy.stats as stats
from tabulate import tabulate
from scipy.stats import chi2, ksone


def menu():
    cad = '\nMenu de Opciones\n' \
          '==============================================\n' \
          '1 ----- Distribución Uniforme \n' \
          '2 ----- Distribución Exponencial\n' \
          '3 ----- Distribución Normal\n' \
          '0 ----- Salir\n'
    return print(cad)


def generar_uniforme(a, b, tam_muestra):
    # Generar variables aleatorias con distribucion uniforme utilzando la formula 𝑋 = 𝐴 + 𝑅𝑁𝐷(𝐵 − 𝐴)
    list_rand = []
    for i in range(tam_muestra):
        ran = random.uniform(0, 1)
        list_rand.append(a + ran * (b - a))
    return list_rand
    # return [round(random.uniform(a, b), 4) for _ in range(tam_muestra)]


def generar_exponencial(lam, tam_muestra):
    # Generar variables aleatorias con distribucion exponencial utilzando la formula X = -1/lam * ln(1 - RND)
    list_rand = []
    for i in range(tam_muestra):
        ran = random.uniform(0, 1)
        list_rand.append(-(1 / lam) * (math.log(1 - ran)))
    return list_rand
    # return [round(random.expovariate(1 / lam), 4) for _ in range(tam_muestra)]


def generar_normal(med, desviacion, tam_muestra):
    # Generar variables aleatorias con distribucion normal utilzando el metodo de Box-Muller
    list_rand = []
    for i in range(tam_muestra):
        ran1 = random.uniform(0, 1)
        ran2 = random.uniform(0, 1)
        list_rand.append(((math.sqrt(-2 * math.log(ran1))) * math.cos(2 * math.pi * ran2)) * desviacion + med)
    return list_rand


def elegirIntervalo():
    num = int(input("Ingrese el número de intervalos (10, 12, 16 o 23): "))
    while num not in (10, 12, 16, 23):
        print("Cantidad de intervalos incorrecta")
        num = int(input("Ingrese el número de intervalos (10, 12, 16 o 23): "))
    return num


def calcular_limites(minimo, maximo, amplitud, intervalos):
    # Calcular los limites inferiores (minimo + amplitud) de todos los intervalos
    limites_inferiores = [minimo + i * amplitud for i in range(intervalos)]

    # Calcular los limites superiores (limite_inferior + amplitud) de todos los intervalos
    limites_superiores = [limite_inf + amplitud for limite_inf in limites_inferiores]

    # Modifico el ultimo limite superior
    limites_superiores[-1] = maximo
    return limites_inferiores, limites_superiores


def calcular_frec_obs(datos, minimo, amplitud, intervalos):
    frecuencias = [0] * intervalos

    for dato in datos:
        intervalo = int((dato - minimo) / amplitud)
        if intervalo == intervalos:
            intervalo -= 1  # Asegurarse de que los valores que caen en el límite superior se incluyan en el último
            # intervalo
        frecuencias[intervalo] += 1
    return frecuencias


def calcular_frec_esp_uniforme(n, intervalos):
    # La frecuencia esperada es la misma para todos los intervalos
    return [n / intervalos] * intervalos


def calcular_frec_esp_exponencial(intervalos, lim_inf, lim_sup, lam, n):
    # La frecuencia esperada es mayor en los primeros intervalos y disminuye en los siguientes intervalos
    frecuencia_esp = []

    for i in range(intervalos):
        frec_esp = ((1 - math.exp(-lam * lim_sup[i])) - (1 - math.exp(-lam * lim_inf[i]))) * n
        frecuencia_esp.append(frec_esp)

    return frecuencia_esp


def calcular_frec_esp_normal(intervalos, lim_inf, lim_sup, med, desviacion, n):
    # La frecuencia esperada es mayor en la media y disminuye para intervalos anteriores y posteriores a la media
    frecuencia_esp = []

    for i in range(intervalos):
        cal_sup = stats.norm.cdf(lim_sup[i], loc=med, scale=desviacion)
        cal_inf = stats.norm.cdf(lim_inf[i], loc=med, scale=desviacion)

        frecuencia_esp.append((cal_sup - cal_inf) * n)
    return frecuencia_esp


def calcular_chi_cuadrado(observado, esperado, intervalos):
    # Chi cuadrado va a ser mayor mientras haya mas diferencia entre la frecuencia observada y la frecuencia esperada
    chi = []
    for i in range(intervalos):
        chi.append(((observado[i] - esperado[i]) ** 2) / esperado[i])
    return chi


def calcular_ks(fo, fe, n, intervalos):
    prob_fo = []
    prob_fe = []
    prob_fo_ac = []
    prob_fe_ac = []
    ks = []

    for i in range(intervalos):
        prob_fo.append(fo[i] / n)
        prob_fe.append(fe[i] / n)
        prob_fo_ac.append(sum(fo[0:i + 1]) / n)
        prob_fe_ac.append(sum(fe[0:i + 1]) / n)
        ks.append(abs(prob_fo_ac[i] - prob_fe_ac[i]))

    return prob_fo, prob_fe, prob_fo_ac, prob_fe_ac, ks


def acomodar_frec(matriz, intervalos):
    # Cantidad de filas de la matriz
    num_filas = intervalos

    i = 0

    # Iteracion sobre las filas de la matriz
    while i < num_filas:

        # Verificar si la frecuencia esperada es menor a 5
        if matriz['Frec Esp'][i] < 5:

            # Verificar si es la ultima fila de la tabla
            if i == num_filas - 1:

                # Si el tamaño de la muestra es muy pequeño puede que la frecuencia esperada nunca sea mayor a 5
                if i == 0:
                    break

                # Si es un array, obtener el mínimo y máximo, sino tomar el valor directamente
                if isinstance(matriz['Intervalo'][i - 1], list):
                    int_min = min(matriz['Intervalo'][i - 1])
                else:
                    int_min = matriz['Intervalo'][i - 1]

                if isinstance(matriz['Intervalo'][i], list):
                    int_max = max(matriz['Intervalo'][i])
                else:
                    int_max = matriz['Intervalo'][i]

                # Sumar todos los datos de la tabla con el intervalo anterior en caso de que este en el ultimo intervalo
                # y la frecuencia esperada sea menor a 5

                matriz['Intervalo'][i - 1] = [int_min, int_max]
                matriz['Lim Sup'][i - 1] = matriz['Lim Sup'][i]
                matriz['Frec Obs'][i - 1] += matriz['Frec Obs'][i]
                matriz['Frec Esp'][i - 1] += matriz['Frec Esp'][i]
                matriz['Chi 2'][i - 1] += matriz['Chi 2'][i]

                # Eliminar la fila actual(ultima fila) ya que la sume con la anterior
                del matriz['Intervalo'][i]
                del matriz['Lim Inf'][i]
                del matriz['Lim Sup'][i]
                del matriz['Frec Obs'][i]
                del matriz['Frec Esp'][i]
                del matriz['Chi 2'][i]

                # Salir del while debido a que no tengo mas intervalos
                break

            # La frecuencia esperada es menor a 5 pero no estoy en el ultimo intervalo
            else:
                # Inicializo variable j que es la que me va a determinar cuantos intervalos debo unir
                j = i + 1

                # Busco algun rango de intervalos donde la suma de las frecuencias esperadas sea mayor a 5
                while sum(matriz['Frec Esp'][i:j]) < 5:
                    j += 1

                    # En caso de llegar al final de la tabla y la suma de las frecuencias esperadas es menor que 5
                    # finaliza el while
                    if j == (intervalos - 1):
                        break

                # Si es un array, obtener el mínimo y máximo, sino tomar el valor directamente
                if isinstance(matriz['Intervalo'][i], list):
                    int_min = min(matriz['Intervalo'][i])
                else:
                    int_min = matriz['Intervalo'][i]

                if isinstance(matriz['Intervalo'][i:j][-1], list):
                    int_max = max(matriz['Intervalo'][i:j][-1])
                else:
                    int_max = matriz['Intervalo'][i:j][-1]

                # Sumar todos los datos de la tabla que esten entre i y j
                # El rango del intervalo va a estar definido por el intervalo menor y el intervalo mayor
                matriz['Intervalo'][i] = [int_min, int_max]
                matriz['Lim Inf'][i] = min(matriz['Lim Inf'][i:j])
                matriz['Lim Sup'][i] = max(matriz['Lim Sup'][i:j])
                matriz['Frec Obs'][i] = sum(matriz['Frec Obs'][i:j])
                matriz['Frec Esp'][i] = sum(matriz['Frec Esp'][i:j])
                matriz['Chi 2'][i] = sum(matriz['Chi 2'][i:j])

                # Eliminar las filas que esten entre i y j
                del matriz['Intervalo'][i + 1:j]
                del matriz['Lim Inf'][i + 1:j]
                del matriz['Lim Sup'][i + 1:j]
                del matriz['Frec Obs'][i + 1:j]
                del matriz['Frec Esp'][i + 1:j]
                del matriz['Chi 2'][i + 1:j]

                # Actualizar el número de filas
                num_filas = len(matriz['Frec Esp'])

                # Volver a verificar si la fila actual tiene una frecuencia esperada mayor a 5
                i -= 1

        i += 1

    return matriz


def graficar(datos, intervalos):
    plt.hist(datos, bins=intervalos, edgecolor='blue')
    plt.title('Histograma de Frecuencias')
    plt.xlabel('Intervalos')
    plt.ylabel('Frecuencia')
    plt.grid(False)
    plt.show()


def formatearMatriz(matriz, tieneKS):
    matriz['Lim Inf'] = [f"{x:.4f}" for x in matriz['Lim Inf']]
    matriz['Lim Sup'] = [f"{x:.4f}" for x in matriz['Lim Sup']]
    matriz['Frec Obs'] = [f"{x:.4f}" for x in matriz['Frec Obs']]
    matriz['Frec Esp'] = [f"{x:.4f}" for x in matriz['Frec Esp']]
    matriz['Chi 2'] = [f"{x:.4f}" for x in matriz['Chi 2']]

    if(tieneKS):
        matriz['Prob obs'] = [f"{x:.4f}" for x in matriz['Prob obs']]
        matriz['Prob esp'] = [f"{x:.4f}" for x in matriz['Prob esp']]
        matriz['Prob obs ac'] = [f"{x:.4f}" for x in matriz['Prob obs ac']]
        matriz['Prob esp ac'] = [f"{x:.4f}" for x in matriz['Prob esp ac']]
        matriz['KS'] = [f"{x:.4f}" for x in matriz['KS']]

    return matriz


def main():
    opcion = -1
    datos = []
    xlambda = None
    desviacion = None
    frec_esp = None
    media = None
    nivel_confianza = 0.9

    while opcion != 0:
        menu()

        opcion = int(input("Ingrese su opcion: "))
        if opcion == 0:
            break
        while opcion not in (1, 2, 3):
            print("Opcion incorrecta")
            opcion = int(input("Ingrese su opcion: "))

        n = int(input("Ingrese el tamaño de muestra (MAX 1.000.000): "))
        while n > 1000000:
            print("El tamaño de muestra es muy grande")
            n = int(input("Ingrese el tamaño de muestra (MAX 1.000.000): "))

        intervalos = elegirIntervalo()

        if opcion == 1:
            a = int(input("Ingrese el limite inferior: "))
            b = int(input("Ingrese el limite superior: "))

            while a >= b:
                print("El valor de a debe ser menor que el valor de b")
                a = int(input("Ingrese el a: "))
                b = int(input("Ingrese el b: "))

            datos = generar_uniforme(a, b, n)

        elif opcion == 2:
            xlambda = float(input("Ingrese lambda: "))

            datos = generar_exponencial(xlambda, n)

        elif opcion == 3:
            media = float(input("Ingrese la media: "))
            desviacion = float(input("Ingrese la desviación: "))

            datos = generar_normal(media, desviacion, n)

        # Mostrar datos generados
        print("\nSerie de números generada:")
        print(datos)

        minimo = min(datos)
        maximo = max(datos)
        rango = maximo - minimo
        amplitud = rango / intervalos

        lim_inf, lim_sup = calcular_limites(minimo, maximo, amplitud, intervalos)

        frec_obs = calcular_frec_obs(datos, minimo, amplitud, intervalos)

        if opcion == 1:
            frec_esp = calcular_frec_esp_uniforme(n, intervalos)
        elif opcion == 2:
            frec_esp = calcular_frec_esp_exponencial(intervalos, lim_inf, lim_sup, xlambda, n)
        elif opcion == 3:
            frec_esp = calcular_frec_esp_normal(intervalos, lim_inf, lim_sup, media, desviacion, n)

        chi_cuadrado = calcular_chi_cuadrado(frec_obs, frec_esp, intervalos)

        prob_fo, prob_fe, prob_fo_ac, prob_fe_ac, ks = calcular_ks(frec_obs, frec_esp, n, intervalos)

        # Generar un matriz con todos los datos
        matriz = {'Intervalo': list(range(1, intervalos + 1)),
                  'Lim Inf': lim_inf,
                  'Lim Sup': lim_sup,
                  'Frec Obs': frec_obs,
                  'Frec Esp': frec_esp,
                  'Chi 2': chi_cuadrado,
                  'Prob obs': prob_fo,
                  'Prob esp': prob_fe,
                  'Prob obs ac': prob_fo_ac,
                  'Prob esp ac': prob_fe_ac,
                  'KS': ks
                  }

        matrizChi2 = {'Intervalo': list(range(1, intervalos + 1)),
                      'Lim Inf': lim_inf.copy(),
                      'Lim Sup': lim_sup.copy(),
                      'Frec Obs': frec_obs.copy(),
                      'Frec Esp': frec_esp.copy(),
                      'Chi 2': chi_cuadrado.copy(),
                      }

        matrizAcumulada = acomodar_frec(matrizChi2, intervalos)
        matrizAcumulada['Chi 2'] = calcular_chi_cuadrado(matrizAcumulada['Frec Obs'],
                                                         matrizAcumulada['Frec Esp'],
                                                         len(matrizAcumulada['Frec Esp']))

        chi_tabulado = chi2.ppf(nivel_confianza, intervalos - (opcion - 1) - 1)
        chi2_calc = sum(matrizAcumulada['Chi 2'])

        if (n < 35):
            ks_tabulado = ksone.ppf(nivel_confianza + 0.05, n)
        else:
            ks_tabulado = 1.22 / math.sqrt(n)
        ks_max = max(matriz['KS'])

        matriz = formatearMatriz(matriz, True)
        print(tabulate(matriz, headers="keys", tablefmt="double_grid", numalign="center"))

        matrizAcumulada = formatearMatriz(matrizAcumulada, False)
        print(tabulate(matrizAcumulada, headers="keys", tablefmt="double_grid", numalign="center"))

        if chi2_calc <= chi_tabulado:
            print()
            print("Chi calculado < Chi tabulado")
            print(round(chi2_calc, 4), "<", round(chi_tabulado, 4))
            print("Se acepta la hipótesis nula H0")

        else:
            print()
            print("Chi calculado > Chi tabulado")
            print(round(chi2_calc, 4), ">", round(chi_tabulado, 4))
            print("No se acepta la hipótesis nula H0")

        if ks_max <= ks_tabulado:
            print()
            print("KS calculado < KS tabulado")
            print(round(ks_max, 4), "<", round(ks_tabulado, 4))
            print("Se acepta la hipótesis nula H0")
        else:
            print()
            print("KS calculado > KS tabulado")
            print(round(ks_max, 4), ">", round(ks_tabulado, 4))
            print("No se acepta la hipótesis nula H0")

        graficar(datos, intervalos)


if __name__ == "__main__":
    main()
