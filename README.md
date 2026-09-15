# TP Integrador — Módulo R
Maestría en Econometría — Universidad Torcuato Di Tella
Laboratorio de Programación en Python y R — Docente: Ian Evangelos Bounos

## Dataset

**Encuesta Permanente de Hogares (EPH)**, INDEC — microdatos individuales, **3er trimestre de 2025**, aglomerados urbanos.
Fuente: https://www.indec.gob.ar/indec/web/Institucional-Indec-BasesDeDatos

Archivos usados (incluidos en `data/`):
- `personas_tot_urb_3T_2025.txt` — base individual (71.322 personas, 204 variables), usada como fuente principal.
- `hogares_tot_urb_3T_2025.txt` — base de hogares (25.190 hogares), no se usa en el modelo final pero se dejó disponible para variables a nivel hogar.

## Pregunta y técnica

**Pregunta:** entre las personas ocupadas, ¿qué factores explican el ingreso de la ocupación principal? (edad/experiencia, sexo, nivel educativo, categoría ocupacional, jornada laboral).

**Técnica:** regresión lineal sobre `log(ingreso)` con una ecuación de Mincer extendida, complementada con:
- selección de variables por AIC/BIC (`step()`)
- regularización Ridge y LASSO (`glmnet`)

Se eligió esta técnica (y no series de tiempo o datos de panel) porque la EPH de un solo trimestre es un corte transversal, sin dimensión temporal repetida.

Muestra final: personas ocupadas, entre 18 y 65 años, con ingreso de la ocupación principal válido y positivo (se excluyen los códigos `-9` de "no sabe/no contesta" en ingreso y `-1` en edad, y `999` en horas trabajadas). Detalle completo de filtros y decisiones sobre valores atípicos en el script y en el informe.

## Cómo correr el proyecto

1. Clonar el repositorio.
2. Abrir `scripts/tp_integrador_R.R` en RStudio (o correr con `Rscript`). El script usa rutas relativas: debe ejecutarse desde la raíz del repositorio, con los archivos de `data/` en su lugar.
3. Instalar dependencias si hace falta (primera vez):
   ```r
   install.packages(c("tidyverse", "broom", "lmtest", "sandwich", "car", "glmnet"))
   ```
4. Correr el script de punta a punta. Usa `set.seed(2026)` para reproducibilidad en los splits train/test y en la cross-validation.

El script está organizado en tres partes que corresponden al enunciado del TP:
- **Parte 1** — carga de datos y justificación de dataset/técnica.
- **Parte 2** — análisis exploratorio (EDA).
- **Parte 3** — modelo, diagnóstico, validación y regularización.

## Estructura del repositorio

```
├── README.md
├── data/
│   ├── personas_tot_urb_3T_2025.txt
│   └── hogares_tot_urb_3T_2025.txt
├── scripts/
│   └── tp_integrador_R.R
└── informe/
    └── informe_TP_R.pdf
```

## Informe

El informe completo (sin código, con gráficos e interpretación) está en `informe/Informe_TP_R.pdf`.
