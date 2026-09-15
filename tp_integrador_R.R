# =============================================================================
# Trabajo Practico Integrador - Modulo R
# Maestria en Econometria - Universidad Torcuato Di Tella
# Docente: Ian Evangelos Bounos | Ano 2026
#
# Tema: Determinantes del ingreso laboral
# Dataset: EPH (INDEC) - microdatos individuales, 3er trimestre 2025
# Tecnica: Regresion lineal (ecuacion de Mincer extendida) + seleccion de
#          variables (step) + regularizacion (Ridge / LASSO)
#
# Cada bloque esta identificado con la Parte del enunciado que responde
# (Parte 1: dataset y tecnica | Parte 2: EDA | Parte 3: modelo e interpretacion)
# =============================================================================

set.seed(2026)

library(tidyverse)
library(broom)
library(lmtest)
library(sandwich)
library(car)      # vif(), linearHypothesis()
library(glmnet)   # Ridge y LASSO

# =============================================================================
# PARTE 1 - Eleccion del dataset y de la tecnica
# =============================================================================
#
# Dataset elegido: Encuesta Permanente de Hogares (EPH, INDEC), base individual,
# 3er trimestre de 2025. Fuente: https://www.indec.gob.ar (microdatos EPH).
#
# Pregunta que nos interesa responder: entre las personas ocupadas, que
# factores (edad/experiencia, sexo, nivel educativo, categoria ocupacional,
# jornada laboral) explican el ingreso de la ocupacion principal?
#
# Tecnica elegida: regresion lineal sobre log(ingreso) siguiendo una ecuacion
# de Mincer extendida, con seleccion de variables por AIC/BIC (step()) y
# regularizacion (Ridge y LASSO via glmnet) para evaluar si conviene
# simplificar el modelo. Es coherente con el dataset porque la variable
# dependiente es continua (no es una serie de tiempo ni un panel: la EPH de
# un solo trimestre es un corte transversal).
#
# ---------------------------------------------------------------------------
# 1.1 - Carga de datos
# ---------------------------------------------------------------------------

personas <- read_delim(
  "personas_tot.urb_3T_2025.txt",
  delim = ";",
  locale = locale(decimal_mark = "."),
  show_col_types = FALSE
)

cat("Personas cargadas:", nrow(personas), "observaciones,", ncol(personas), "variables\n")

# =============================================================================
# PARTE 2 - Analisis exploratorio de datos (EDA)
# =============================================================================
#
# 2.1 - Estructura general
# ---------------------------------------------------------------------------

cat("\n", strrep("=", 70), "\n")
cat("PARTE 2 - EDA: estructura general\n")
cat(strrep("=", 70), "\n")

cat("Dimensiones:", nrow(personas), "filas x", ncol(personas), "columnas\n")
cat("Cobertura: EPH 3er trimestre 2025, aglomerados urbanos\n")

# Distribucion de la variable ESTADO (condicion de actividad)
# 1 = ocupado, 2 = desocupado, 3 = inactivo, 4 = menor de 10 anos, 0 = no entrevistado
personas |>
  count(ESTADO) |>
  mutate(pct = round(100 * n / sum(n), 1))

# ---------------------------------------------------------------------------
# 2.2 - Valores faltantes y codigos especiales (previo a cualquier filtro)
# ---------------------------------------------------------------------------
#
# La EPH no usa NA sino codigos numericos para "no sabe/no contesta" u otras
# situaciones. Antes de tratarlos como faltantes hay que identificarlos:
#   - P21 == -9        -> ingreso "no sabe/no contesta" (no es un ingreso real)
#   - CH06 == -1       -> edad no valida
#   - PP3E_TOT == 999  -> horas trabajadas "no sabe/no contesta"
#   - PP04D_COD        -> ~56% de faltantes reales (codigo de actividad),
#                          no se usa en el modelo por su alta tasa de faltantes

faltantes <- tibble(
  variable = c("P21 (ingreso ppal.)", "CH06 (edad)", "PP3E_TOT (horas)"),
  codigo_especial = c(-9, -1, 999),
  n_casos = c(
    sum(personas$P21 == -9),
    sum(personas$CH06 == -1),
    sum(personas$PP3E_TOT == 999)
  )
) |>
  mutate(pct = round(100 * n_casos / nrow(personas), 2))

print(faltantes)

# ---------------------------------------------------------------------------
# 2.3 - Filtro de la muestra y construccion de variables
# ---------------------------------------------------------------------------
#
# Nos quedamos con personas ocupadas (ESTADO == 1), en edad activa (18-65),
# con ingreso de la ocupacion principal valido y positivo. Se excluyen los
# codigos -9 (P21) y -1 (CH06) por ser "no sabe/no contesta", no ingresos ni
# edades reales. PP3E_TOT == 999 se recodifica como NA por el mismo motivo.
# El nivel educativo "sin instruccion" (NIVEL_ED == 7) se agrupa junto con
# primario, porque tiene muy pocos casos e ingresos similares a ese grupo.

eph <- personas |>
  filter(ESTADO == 1, P21 > 0, CH06 >= 18, CH06 <= 65) |>
  mutate(
    log_ingreso = log(P21),
    sexo = factor(CH04, levels = c(1, 2), labels = c("Varon", "Mujer")),
    edad = CH06,
    experiencia = pmax(edad - 18, 0),
    exper2 = experiencia^2,
    nivel_ed = case_when(
      NIVEL_ED %in% c(1, 2, 7) ~ "Primario o menos",
      NIVEL_ED %in% c(3, 4)    ~ "Secundario",
      NIVEL_ED %in% c(5, 6)    ~ "Superior",
      TRUE ~ NA_character_
    ),
    nivel_ed = factor(nivel_ed,
                       levels = c("Primario o menos", "Secundario", "Superior")),
    cat_ocup = case_when(
      CAT_OCUP == 3 ~ "Obrero/Empleado",
      CAT_OCUP == 1 ~ "Patron",
      CAT_OCUP == 2 ~ "Cuentapropista",
      TRUE ~ NA_character_
    ),
    cat_ocup = factor(cat_ocup,
                       levels = c("Obrero/Empleado", "Patron", "Cuentapropista")),
    horas_trab = ifelse(PP3E_TOT >= 999, NA, PP3E_TOT),
    jornada_comp = factor(ifelse(horas_trab >= 35, "Completa", "Parcial"))
  ) |>
  filter(!is.na(nivel_ed), !is.na(cat_ocup), !is.na(horas_trab))

cat("\nMuestra final para el modelo:", nrow(eph), "observaciones\n")
cat("(", nrow(personas), "personas totales ->", sum(personas$ESTADO == 1),
    "ocupadas ->", nrow(eph), "con datos completos y validos )\n")

# ---------------------------------------------------------------------------
# 2.4 - Distribuciones univariadas
# ---------------------------------------------------------------------------

# Distribucion del ingreso (nivel vs log) - justifica trabajar en logaritmos
p_ingreso_nivel <- eph |>
  ggplot(aes(x = P21)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.8) +
  labs(title = "Distribucion del ingreso de la ocupacion principal",
       subtitle = "Fuertemente sesgada a la derecha",
       x = "Ingreso ($)", y = "Frecuencia") +
  theme_minimal(base_size = 13)

p_ingreso_log <- eph |>
  ggplot(aes(x = log_ingreso)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.8) +
  labs(title = "Distribucion del log(ingreso)",
       subtitle = "Mucho mas cercana a una normal -> justifica log(P21) como variable dependiente",
       x = "log(Ingreso)", y = "Frecuencia") +
  theme_minimal(base_size = 13)

print(p_ingreso_nivel)
print(p_ingreso_log)

# Edad / experiencia
p_edad <- eph |>
  ggplot(aes(x = edad)) +
  geom_histogram(bins = 40, fill = "darkorange", alpha = 0.8) +
  labs(title = "Distribucion de la edad (muestra ocupados 18-65)",
       x = "Edad", y = "Frecuencia") +
  theme_minimal(base_size = 13)

print(p_edad)

# Composicion de la muestra por sexo, nivel educativo y categoria ocupacional
p_composicion <- eph |>
  count(sexo, nivel_ed, cat_ocup) |>
  ggplot(aes(x = nivel_ed, y = n, fill = sexo)) +
  geom_col(position = "dodge") +
  facet_wrap(~cat_ocup) +
  labs(title = "Composicion de la muestra",
       x = "Nivel educativo", y = "Cantidad de personas", fill = "Sexo") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

print(p_composicion)

# ---------------------------------------------------------------------------
# 2.5 - Relaciones bivariadas relevantes para el modelo
# ---------------------------------------------------------------------------

# Ingreso segun sexo
p_sexo <- eph |>
  ggplot(aes(x = sexo, y = log_ingreso, fill = sexo)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "log(Ingreso) segun sexo",
       subtitle = "Brecha visible a favor de los varones",
       x = NULL, y = "log(Ingreso)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

print(p_sexo)

# Ingreso segun nivel educativo
p_niveled <- eph |>
  ggplot(aes(x = nivel_ed, y = log_ingreso, fill = nivel_ed)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "log(Ingreso) segun nivel educativo",
       x = NULL, y = "log(Ingreso)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

print(p_niveled)

# Ingreso segun categoria ocupacional
p_catocup <- eph |>
  ggplot(aes(x = cat_ocup, y = log_ingreso, fill = cat_ocup)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "log(Ingreso) segun categoria ocupacional",
       x = NULL, y = "log(Ingreso)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

print(p_catocup)

# Ingreso vs experiencia (perfil edad-ingreso, forma de U invertida esperada)
p_experiencia <- eph |>
  ggplot(aes(x = experiencia, y = log_ingreso)) +
  geom_point(alpha = 0.15, size = 0.8, color = "steelblue") +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(title = "log(Ingreso) vs anos de experiencia potencial",
       subtitle = "La curva sugiere una relacion no lineal (justifica incluir experiencia^2)",
       x = "Experiencia potencial (anos)", y = "log(Ingreso)") +
  theme_minimal(base_size = 13)

print(p_experiencia)

# ---------------------------------------------------------------------------
# 2.6 - Valores atipicos: que se hizo y por que
# ---------------------------------------------------------------------------
#
# Se identificaron y trataron los siguientes casos atipicos / codigos
# especiales (a documentar en el informe, Parte 2):
#   - P21 == -9 (no sabe/no contesta): EXCLUIDO, no es un ingreso real.
#   - CH06 == -1 (edad invalida): no aparece en la submuestra de ocupados
#     18-65, pero se documenta que existe en la base completa.
#   - PP3E_TOT == 999 (horas ns/nc, 8 casos): recodificado a NA y excluido.
#   - PP3E_TOT == 0 (587 casos en la base de ocupados): personas con licencia
#     o vacaciones esa semana, que igual declaran ingreso. Se DEJAN en la
#     muestra (no es un error, es una situacion real), pero se marca como
#     jornada "Parcial" y se comenta la limitacion en el informe.
#   - Ingresos extremadamente altos (P21 varios ordenes de magnitud por
#     encima de la mediana): no se recortan a mano; usar log(ingreso) ya
#     reduce mucho su influencia, y se vuelve a revisar en el diagnostico
#     del modelo (Parte 3, VIF / distancia de Cook).

cat("\nCasos con PP3E_TOT == 0 en la muestra final:",
    sum(eph$horas_trab == 0), "\n")

# =============================================================================
# PARTE 3 - Aplicacion de la tecnica e interpretacion
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("PARTE 3 - Modelo de determinantes del ingreso\n")
cat(strrep("=", 70), "\n")

# ---------------------------------------------------------------------------
# 3.1 - Modelos anidados (ecuacion de Mincer extendida)
# ---------------------------------------------------------------------------

m0 <- lm(log_ingreso ~ experiencia + exper2 + nivel_ed, data = eph)
m1 <- lm(log_ingreso ~ experiencia + exper2 + nivel_ed + sexo, data = eph)
m2 <- lm(log_ingreso ~ experiencia + exper2 + nivel_ed + sexo +
           cat_ocup + jornada_comp, data = eph)

summary(m2)

# Test F: agregar sexo, y luego categoria ocupacional + jornada, mejora el ajuste?
anova(m0, m1)
anova(m1, m2)

# Hipotesis especifica: existe brecha salarial de genero?
linearHypothesis(m2, "sexoMujer = 0")
# Con errores robustos HC3 (recomendado si hay heterocedasticidad, ver 3.2)
linearHypothesis(m2, "sexoMujer = 0", vcov = vcovHC(m2, type = "HC3"))

# ---------------------------------------------------------------------------
# 3.2 - Diagnostico del modelo
# ---------------------------------------------------------------------------

par(mfrow = c(2, 2))
plot(m2)
par(mfrow = c(1, 1))

vif(m2)
bptest(m2)

# Si bptest() rechaza homocedasticidad, reportar con errores robustos:
coeftest(m2, vcov = vcovHC(m2, type = "HC3")) |>
  tidy() |>
  mutate(across(where(is.numeric), ~round(., 4)))

# ---------------------------------------------------------------------------
# 3.3 - Validacion: train/test split
# ---------------------------------------------------------------------------

n <- nrow(eph)
idx_train <- sample(1:n, size = round(0.7 * n), replace = FALSE)
idx_test <- setdiff(1:n, idx_train)
eph_train <- eph[idx_train, ]
eph_test <- eph[idx_test, ]

m_train <- lm(log_ingreso ~ experiencia + exper2 + nivel_ed + sexo +
                cat_ocup + jornada_comp, data = eph_train)
pred_test <- predict(m_train, newdata = eph_test)

ss_res <- sum((eph_test$log_ingreso - pred_test)^2)
ss_tot <- sum((eph_test$log_ingreso - mean(eph_test$log_ingreso))^2)
r2_test <- 1 - ss_res / ss_tot
rmse_test <- sqrt(mean((eph_test$log_ingreso - pred_test)^2))

cat(sprintf("\nR2 train: %.4f\n", summary(m_train)$r.squared))
cat(sprintf("R2 test:  %.4f\n", r2_test))
cat(sprintf("RMSE test: %.4f\n", rmse_test))

p_pred_obs <- data.frame(observado = eph_test$log_ingreso, predicho = pred_test) |>
  ggplot(aes(x = predicho, y = observado)) +
  geom_point(alpha = 0.2, size = 0.8, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Predichos vs observados - conjunto de test",
       subtitle = sprintf("R2 test = %.3f | RMSE = %.3f", r2_test, rmse_test),
       x = "log(Ingreso) predicho", y = "log(Ingreso) observado") +
  theme_minimal(base_size = 13)

print(p_pred_obs)

# ---------------------------------------------------------------------------
# 3.4 - Seleccion de variables por AIC/BIC
# ---------------------------------------------------------------------------

m_min <- lm(log_ingreso ~ 1, data = eph)
m_max <- lm(log_ingreso ~ experiencia + exper2 + nivel_ed + sexo +
              cat_ocup + jornada_comp, data = eph)

m_stepwise <- step(m_min, scope = list(lower = m_min, upper = m_max),
                    direction = "both", trace = 0)
cat("\nVariables seleccionadas (stepwise AIC):\n")
print(formula(m_stepwise))

m_bic <- step(m_min, scope = list(lower = m_min, upper = m_max),
              direction = "both", k = log(nrow(eph)), trace = 0)
cat("\nVariables seleccionadas (stepwise BIC):\n")
print(formula(m_bic))

# ---------------------------------------------------------------------------
# 3.5 - Regularizacion: Ridge y LASSO (glmnet)
# ---------------------------------------------------------------------------
#
# Con pocas variables (6 predictores, todas sustantivas) no esperamos que
# Ridge/LASSO cambien mucho el resultado frente a OLS; el objetivo aca es
# comparar y mostrar que el modelo es estable (o, si LASSO descarta alguna
# categoria, discutir por que).

X <- model.matrix(log_ingreso ~ experiencia + exper2 + nivel_ed + sexo +
                     cat_ocup + jornada_comp, data = eph)[, -1]
y <- eph$log_ingreso

X_tr <- X[idx_train, ]; y_tr <- y[idx_train]
X_te <- X[idx_test, ];  y_te <- y[idx_test]

ridge_cv <- cv.glmnet(X_tr, y_tr, alpha = 0, nfolds = 10)
lasso_cv <- cv.glmnet(X_tr, y_tr, alpha = 1, nfolds = 10)

pred_ridge <- predict(ridge_cv, s = "lambda.min", newx = X_te)
pred_lasso <- predict(lasso_cv, s = "lambda.min", newx = X_te)

rmse_ridge <- sqrt(mean((y_te - pred_ridge)^2))
rmse_lasso <- sqrt(mean((y_te - pred_lasso)^2))

comparacion_rmse <- tibble(
  Metodo = c("OLS", "Ridge (lambda.min)", "LASSO (lambda.min)"),
  RMSE = c(rmse_test, rmse_ridge, rmse_lasso)
) |>
  mutate(RMSE = round(RMSE, 4)) |>
  arrange(RMSE)

print(comparacion_rmse)

# Variables que LASSO descarta (coeficiente == 0)
coef_lasso <- coef(lasso_cv, s = "lambda.min")
vars_fuera <- rownames(coef_lasso)[coef_lasso[, 1] == 0]
cat("\nVariables descartadas por LASSO (coeficiente = 0):",
    ifelse(length(vars_fuera) == 0, "ninguna", paste(vars_fuera, collapse = ", ")),
    "\n")

# Comparacion visual de coeficientes OLS vs Ridge vs LASSO
coef_ols <- coef(m2)[-1]
coef_ridge_v <- as.vector(coef(ridge_cv, s = "lambda.min"))[-1]
coef_lasso_v <- as.vector(coef(lasso_cv, s = "lambda.min"))[-1]

comp_coefs <- tibble(
  variable = names(coef_ols),
  OLS = coef_ols,
  Ridge = coef_ridge_v,
  LASSO = coef_lasso_v
) |>
  pivot_longer(-variable, names_to = "metodo", values_to = "coef") |>
  mutate(metodo = factor(metodo, levels = c("OLS", "Ridge", "LASSO")))

p_comp_coefs <- comp_coefs |>
  ggplot(aes(x = reorder(variable, abs(coef)), y = coef, fill = metodo)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0, color = "gray50") +
  coord_flip() +
  labs(title = "Comparacion de coeficientes: OLS vs Ridge vs LASSO",
       x = NULL, y = "Coeficiente estimado", fill = "Metodo") +
  theme_minimal(base_size = 12)

print(p_comp_coefs)

# ---------------------------------------------------------------------------
# 3.6 - Interpretacion (guia para el informe en PDF - Parte 3)
# ---------------------------------------------------------------------------
#
# En el informe conviene cubrir, con los resultados de arriba:
#   - Magnitud e interpretacion de los coeficientes de m2 (en % aprox., ya
#     que la variable dependiente esta en logaritmos): sexo, nivel educativo,
#     categoria ocupacional, jornada.
#   - Resultado del test de brecha de genero (linearHypothesis) y si sigue
#     siendo significativo con errores robustos.
#   - Que dice el diagnostico (VIF, Breusch-Pagan): hay heterocedasticidad?
#     hay multicolinealidad relevante (aparte de experiencia/experiencia^2)?
#   - Que tan bien generaliza el modelo (R2/RMSE en test).
#   - Si Ridge/LASSO cambian la conclusion o confirman el modelo OLS.
#   - Limitaciones: es un corte transversal (no permite hablar de causalidad
#     ni de cambios en el tiempo), variables omitidas (rama de actividad,
#     region, antiguedad real en el empleo en vez de experiencia potencial).
