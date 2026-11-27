"""Calculadora eTFG — educational implementation of published CKD-EPI equations."""
import math
from decimal import Decimal, ROUND_HALF_UP

VERSION = '1.0'
UPDATED = '26/11/2026'
AUTHORS = ['Abel Suárez-Castro', 'Carmen I. Bárcenas-López', 'Fernando Rosales-López', 'Genaro F. Sánchez-Mejorada']
UNITS = 'mL/min/1.73 m²'
NIDDK = 'https://www.niddk.nih.gov/research-funding/research-programs/kidney-clinical-research-epidemiology/laboratory/glomerular-filtration-rate-equations/adults'
RISK = ((0,1,2),(0,1,2),(1,2,3),(2,3,3),(3,3,3),(3,3,3))
G_NAMES = ('G1','G2','G3a','G3b','G4','G5')
G_RANGES = ('≥90','60–89','45–59','30–44','15–29','<15')
COLORS = ('#dcebd8','#fae7a0','#f7c38b','#efb0ad')


def parse_value(text, label, optional=False, zero=False):
    text = str(text).strip()
    if not text and optional:
        return None
    try:
        value = float(text.replace(',', '.'))
    except ValueError:
        raise ValueError(label) from None
    if not math.isfinite(value) or value < 0 or (value == 0 and not zero):
        raise ValueError(label)
    return value


def reported(value):
    return int(Decimal(str(value)).quantize(Decimal('1'), rounding=ROUND_HALF_UP))


def g_index(value):
    return next((i for i, cutoff in enumerate((90,60,45,30,15)) if value >= cutoff), 5)


def a_index(value, unit):
    if value is None:
        return None
    low, high = (30,300) if unit == 'mg/g' else (3,30)
    return 0 if value < low else 1 if value <= high else 2


def estimate(age, sex, cr=None, cys=None):
    """Return raw results and intermediate terms; cr mg/dL, cys mg/L."""
    if not math.isfinite(age) or age < 18 or sex not in ('F','M'):
        raise ValueError('Adult age and sex required')
    if cr is None and cys is None:
        raise ValueError('At least one marker is required')
    if any(not math.isfinite(v) or v <= 0 for v in (cr,cys) if v is not None):
        raise ValueError('Invalid biomarker')
    female = sex == 'F'
    k = 0.7 if female else 0.9
    results = {}
    def compute(name, base, alpha, cr_exp, cy_low, cy_high, age_base, sf, use_cr, use_cys):
        terms = [('μ',base),('Age / Edad',age)]
        factors = [base]
        if use_cr:
            x = cr/k
            mn,mx = min(x,1),max(x,1)
            p,q = mn**alpha,mx**cr_exp
            terms += [('κ',k),('α',alpha),('SCr (mg/dL)',cr),('SCr / κ',x),('min(SCr / κ, 1)',mn),('max(SCr / κ, 1)',mx),('min(SCr / κ, 1)^α',p),('max(SCr / κ, 1)^'+str(cr_exp),q)]
            factors += [p,q]
        if use_cys:
            y = cys/0.8
            mn,mx = min(y,1),max(y,1)
            p,q = mn**cy_low,mx**cy_high
            terms += [('SCys (mg/L)',cys),('SCys / 0.8',y),('min(SCys / 0.8, 1)',mn),('max(SCys / 0.8, 1)',mx),('min(SCys / 0.8, 1)^'+str(cy_low),p),('max(SCys / 0.8, 1)^'+str(cy_high),q)]
            factors += [p,q]
        age_factor = age_base**age
        factors += [age_factor,sf]
        value = math.prod(factors)
        if not math.isfinite(value) or value <= 0:
            raise ValueError('Result outside numerical range')
        terms += [(str(age_base)+'^Age',age_factor),('F',sf),('eGFR / eTFG',value)]
        results[name] = {'value':value,'terms':terms,'factors':factors}
    if cr is not None:
        compute('cr',142,-.241 if female else -.302,-1.2,None,None,.9938,1.012 if female else 1,True,False)
    if cys is not None:
        compute('cys',133,None,None,-.499,-1.328,.996,.932 if female else 1,False,True)
    if cr is not None and cys is not None:
        compute('both',135,-.219 if female else -.144,-.544,-.323,-.778,.9961,.963 if female else 1,True,True)
    return results


FORMULAS = {
'cr': r'142\,\min(SCr/\kappa,1)^\alpha\,\max(SCr/\kappa,1)^{-1.200}\,0.9938^{Age}\,F',
'cys': r'133\,\min(SCys/0.8,1)^{-0.499}\,\max(SCys/0.8,1)^{-1.328}\,0.996^{Age}\,F',
'both': r'135\,\min(SCr/\kappa,1)^\alpha\,\max(SCr/\kappa,1)^{-0.544}\,\min(SCys/0.8,1)^{-0.323}\,\max(SCys/0.8,1)^{-0.778}\,0.9961^{Age}\,F'}
CONSTANTS = [('μ','142','133','135'),('κ F','0.7','—','0.7'),('κ M','0.9','—','0.9'),('α F','−0.241','—','−0.219'),('α M','−0.302','—','−0.144'),('SCr high / superior','−1.200','—','−0.544'),('SCys knot / inflexión','—','0.8','0.8'),('SCys low / inferior','—','−0.499','−0.323'),('SCys high / superior','—','−1.328','−0.778'),('Age base / edad','0.9938','0.996','0.9961'),('F (female / mujer)','1.012','0.932','0.963'),('F (male / hombre)','1','1','1')]


def main():
    import streamlit as st
    import pandas as pd
    st.set_page_config(page_title='Calculadora eTFG', page_icon='🧪',layout='wide')
    st.markdown('''<style>
    .block-container {max-width:1120px;padding-top:5rem;padding-bottom:2rem;}
    div[data-testid="stMetric"] {border:1px solid #8b9da344;border-radius:12px;padding:16px;}
    .renal-kicker {color:#008b83;font-size:.8rem;letter-spacing:.13em;font-weight:700;}
    .renal-title {font-size:2.5rem;font-weight:750;line-height:1.1;margin:.3rem 0;}
    .renal-sub {font-size:1.05rem;opacity:.82;margin-bottom:1rem;}
    </style>''',unsafe_allow_html=True)
    lang = st.sidebar.radio('Idioma / Language',['Español','English'],key='lang')
    es = lang == 'Español'
    def T(a,b): return a if es else b
    title = T('Calculadora eTFG','eGFR Calculator')
    names = {'cr':T('Creatinina · 2021','Creatinine · 2021'),'cys':T('Cistatina C · 2012','Cystatin C · 2012'),'both':T('Creatinina + cistatina C · 2021','Creatinine + cystatin C · 2021')}
    descriptions = T(['Normal o elevada','Ligeramente disminuida','Ligera a moderadamente disminuida','Moderada a gravemente disminuida','Gravemente disminuida','Falla renal'],['Normal or high','Mildly decreased','Mildly to moderately decreased','Moderately to severely decreased','Severely decreased','Kidney failure'])
    risks = T(['Bajo','Moderado','Alto','Muy alto'],['Low','Moderate','High','Very high'])
    st.sidebar.markdown('**CKD-EPI 2012 / 2021 · KDIGO 2024**')
    st.sidebar.caption(T('Adultos ≥18 años. Sin coeficiente racial.','Adults ≥18 years. No race coefficient.'))
    st.markdown(f'<div class="renal-kicker">{T("MEDICINA · QUÍMICA CLÍNICA","MEDICINE · CLINICAL CHEMISTRY")}</div><div class="renal-title">{title}</div><div class="renal-sub">{T("Tasa de filtración glomerular estimada · Ecuaciones CKD-EPI", "Estimated glomerular filtration rate · CKD-EPI equations")}</div>',unsafe_allow_html=True)
    st.caption('CKD-EPI 2012 / 2021  ·  KDIGO 2024  ·  v'+VERSION)
    st.info(T('Herramienta académica y educativa para adultos. No diagnostica enfermedad renal ni sustituye la valoración médica. No introduzcas nombres, expedientes ni otros datos identificables.', 'Academic and educational tool for adults. It does not diagnose kidney disease or replace medical assessment. Do not enter names, medical record numbers or other identifying information.'))
    tabs = st.tabs(T(['Calculadora','Cálculos intermedios','Constantes y fórmulas','Clasificación KDIGO','Calidad analítica','Referencias'],['Calculator','Intermediate calculations','Constants and equations','KDIGO categories','Analytical quality','References']))
    results = {}; primary = None; gi = ai = None; acr = None; data = None
    def reset():
        for key in ('age','cr','cys','acr'): st.session_state[key]=''
        st.session_state['sex']=''
        st.session_state['cr_unit']='mg/dL'; st.session_state['acr_unit']='mg/g'
        st.session_state['submitted']=None
    def example():
        reset()
        st.session_state.update({'age':'55','sex':'M','cr':'1.5','cys':'1.8','acr':'45','example_requested':True})
    def clear_cr(): st.session_state['cr']=''; st.session_state['submitted']=None
    def clear_acr(): st.session_state['acr']=''; st.session_state['submitted']=None
    with tabs[0]:
        st.subheader(T('Datos de entrada','Input data'), anchor=False)
        left,right = st.columns(2)
        with left:
            age_text=st.text_input(T('Edad (años)','Age (years)'),key='age',placeholder=T('Ingresa la edad','Enter age'))
        with right:
            sex=st.selectbox(T('Sexo utilizado en la ecuación','Sex used by the equation'),['','F','M'],format_func=lambda x: {'':T('Selecciona una opción','Select an option'),'F':T('Femenino','Female'),'M':T('Masculino','Male')}[x],key='sex')
        left,right = st.columns(2)
        with left:
            cr_unit=st.radio(T('Unidad de creatinina','Creatinine unit'),['mg/dL','µmol/L'],horizontal=True,key='cr_unit',on_change=clear_cr)
            cr_text=st.text_input(T('Creatinina sérica','Serum creatinine')+f' ({cr_unit})',key='cr',placeholder=T('Ingresa el valor','Enter value'))
        with right:
            st.caption(T('Introduce al menos un marcador de filtración.','Enter at least one filtration marker.'))
            cys_text=st.text_input(T('Cistatina C sérica (mg/L)','Serum cystatin C (mg/L)'),key='cys',placeholder=T('Ingresa el valor','Enter value'))
        st.caption(T('Para la ecuación combinada, se prefieren resultados de ambos marcadores obtenidos de la misma muestra. Cambiar unidades borra el campo correspondiente.','For the combined equation, both markers should preferably be measured in the same sample. Changing units clears the corresponding field.'))
        with st.expander(T('Albuminuria (opcional)','Albuminuria (optional)'),expanded=True):
            acr_unit=st.radio(T('Unidad de RAC urinaria','Urine ACR unit'),['mg/g','mg/mmol'],horizontal=True,key='acr_unit',on_change=clear_acr)
            acr_text=st.text_input(T('Relación albúmina/creatinina urinaria (RAC)','Urine albumin-to-creatinine ratio (ACR)')+f' ({acr_unit})',key='acr',placeholder=T('Ingresa el valor','Enter value'))
        signature=(age_text,sex,cr_text,cys_text,acr_text,cr_unit,acr_unit)
        b1,b2,b3=st.columns([2,1,1])
        if b1.button(T('Calcular','Calculate'),type='primary',use_container_width=True) or st.session_state.pop('example_requested',False):
            st.session_state['submitted']=signature
        b2.button(T('Ejemplo ficticio','Sample data'),on_click=example,use_container_width=True)
        b3.button(T('Limpiar','Clear'),on_click=reset,use_container_width=True)
        if st.session_state.get('submitted') == signature:
            try:
                age=parse_value(age_text,T('Edad','Age'))
                if age<18: raise ValueError(T('Edad: se requiere ≥18 años','Age: must be ≥18 years'))
                if sex not in ('F','M'): raise ValueError(T('Selecciona el sexo','Select sex'))
                cr=parse_value(cr_text,T('Creatinina','Creatinine'),True)
                cys=parse_value(cys_text,T('Cistatina C','Cystatin C'),True)
                acr=parse_value(acr_text,T('RAC urinaria','Urine ACR'),True,True)
                if cr is None and cys is None: raise ValueError(T('Introduce creatinina o cistatina C','Enter creatinine or cystatin C'))
                cr_mg=cr/(88.4 if cr_unit=='µmol/L' else 1) if cr is not None else None
                results=estimate(age,sex,cr_mg,cys)
                primary='both' if 'both' in results else 'cr' if 'cr' in results else 'cys'
                value=results[primary]['value']; shown=reported(value); gi=g_index(shown); ai=a_index(acr,acr_unit)
                data=(age,sex,cr,cr_mg,cys,acr,cr_unit,acr_unit)
            except (ValueError,OverflowError,ArithmeticError) as err:
                results={}; primary=None; gi=ai=None
                st.error(T('Revisa los datos: ','Check the inputs: ')+str(err))
        elif st.session_state.get('submitted') is not None:
            st.caption(T('Los datos cambiaron. Pulsa Calcular para actualizar.','Inputs changed. Press Calculate to update.'))
        if primary is not None:
            st.divider(); st.subheader(T('Resultados','Results'), anchor=False)
            c1,c2=st.columns([2,1]); c1.metric(names[primary],f'{shown} {UNITS}'); c2.metric(T('Categoría de filtración','GFR category'),G_NAMES[gi])
            st.caption(descriptions[gi]+' · '+T('Categoría basada en el resultado redondeado de laboratorio.','Category based on the rounded laboratory result.'))
            if g_index(value)!=gi:
                st.warning(T('Resultado cercano a un límite: el valor sin redondear corresponde a ', 'Result near a threshold: the unrounded value corresponds to ')+G_NAMES[g_index(value)]+'.')
            if primary=='both': st.caption(T('Se destaca la ecuación combinada cuando ambos marcadores están disponibles (KDIGO 2024, 1.1.2.1).','The combined equation is highlighted when both markers are available (KDIGO 2024, 1.1.2.1).'))
            if shown<60: st.warning(T('eTFG reportada inferior a 60. Requiere valoración clínica; un resultado aislado no establece cronicidad.','Reported eGFR below 60. Clinical assessment is needed; a single result does not establish chronicity.'))
            else: st.info(T('G1 o G2 no excluyen daño renal. La albuminuria y otros hallazgos son relevantes.','G1 or G2 do not exclude kidney damage. Albuminuria and other findings matter.'))
            if shown<15: st.warning(T('Estimación muy baja: busca valoración médica pronta; si es un hallazgo nuevo o hay síntomas graves, atención urgente.','Very low estimate: seek prompt medical assessment; if newly found or severe symptoms are present, seek urgent care.'))
            if 18<=age<=25: st.caption(T('En edades de 18–25 años, NIDDK recomienda también considerar CKiD U25 (no incluida aquí).','For ages 18–25, NIDDK recommends also considering CKiD U25 (not included here).'))
            if len(results)>1:
                st.subheader(T('Comparación de estimaciones','Comparison of estimates'), anchor=False)
                df=pd.DataFrame({T('Ecuación','Equation'):[names[k] for k in results],UNITS:[reported(v['value']) for v in results.values()]})
                st.dataframe(df,hide_index=True,use_container_width=True)
                if 'both' in results:
                    delta=results['cys']['value']-results['cr']['value']
                    st.caption(T('Diferencia cistatina C − creatinina (sin redondear): ','Cystatin C − creatinine difference (unrounded): ')+f'{delta:+.2f} {UNITS} ({100*delta/results["cr"]["value"]:+.1f}% '+T('respecto a creatinina). No es un umbral diagnóstico.','relative to creatinine). Not a diagnostic threshold.'))
            if ai is not None:
                st.metric(T('Categoría de albuminuria','Albuminuria category'),f'A{ai+1}')
                st.caption(T('Consulta la matriz en la pestaña Clasificación KDIGO.','See the matrix in the KDIGO categories tab.'))
    with tabs[1]:
        st.header(T('Procedimiento matemático','Mathematical procedure'), anchor=False)
        st.caption(T('Se conserva precisión interna; los decimales visibles son solo para auditoría.','Internal precision is retained; displayed decimals are for auditing only.'))
        if not results: st.info(T('Calcula primero un ejemplo o introduce tus valores.','Calculate an example or enter your values first.'))
        else:
            if data[2] is not None and data[6]=='µmol/L': st.write(f'SCr = {data[2]:g} / 88.4 = {data[3]:.10g} mg/dL')
            for k,v in results.items():
                with st.expander(names[k],expanded=(k==primary)):
                    st.latex(r'eGFR = '+FORMULAS[k]); st.dataframe(pd.DataFrame({T('Término','Term'):[x[0] for x in v['terms']],T('Valor','Value'):[f'{x[1]:.10g}' for x in v['terms']]}),hide_index=True,use_container_width=True)
                    st.latex(r'\times '.join(f'{x:.10g}' for x in v['factors'])+r'\approx '+f'{v["value"]:.8g}')
                    st.write(T('Resultado reportado: ','Reported result: ')+f'**{reported(v["value"])} {UNITS}**')
    with tabs[2]:
        st.header(T('Ecuaciones y constantes','Equations and constants'), anchor=False)
        st.caption(T('SCr: creatinina sérica en mg/dL. SCys: cistatina C sérica en mg/L. Age: edad en años. F: multiplicador por sexo.','SCr: serum creatinine in mg/dL. SCys: serum cystatin C in mg/L. Age: age in years. F: sex multiplier.'))
        for k in ('cr','cys','both'):
            st.markdown('**'+names[k]+'**'); st.latex('eGFR = '+FORMULAS[k])
        st.dataframe(pd.DataFrame(CONSTANTS,columns=[T('Constante','Constant'),names['cr'],names['cys'],names['both']]),hide_index=True,use_container_width=True)
        st.write(T('min selecciona el menor de dos valores; max selecciona el mayor. κ y 0.8 son puntos de inflexión matemáticos, no límites de normalidad.','min selects the smaller value; max selects the larger. κ and 0.8 are mathematical knots, not normal reference limits.'))
        st.caption(T('Fuentes: Inker 2012, tabla 2; Inker 2021, tabla 2; Miller 2022, tablas 3–4; NIDDK. Para cistatina C 2012 se usa 0.996, tal como en el artículo original y NIDDK.','Sources: Inker 2012, Table 2; Inker 2021, Table 2; Miller 2022, Tables 3–4; NIDDK. For 2012 cystatin C we use 0.996 as in the original paper and NIDDK.'))
    with tabs[3]:
        st.header(T('Clasificación KDIGO 2024','KDIGO 2024 classification'), anchor=False)
        st.write(T('La ERC se define por anomalías de estructura o función renal durante al menos tres meses, con implicaciones para la salud. Se clasifica por causa, categoría G y categoría A (CGA).','CKD is defined by abnormalities of kidney structure or function for at least three months, with implications for health. It is classified by cause, GFR category and albuminuria category (CGA).'))
        st.dataframe(pd.DataFrame({'G':G_NAMES,UNITS:G_RANGES,T('Descripción','Description'):descriptions}),hide_index=True,use_container_width=True)
        st.dataframe(pd.DataFrame({'A':['A1','A2','A3'],'mg/g':['<30','30–300','>300'],'mg/mmol':['<3','3–30','>30'],T('Descripción','Description'):T(['Normal a ligeramente aumentada','Moderadamente aumentada','Gravemente aumentada'],['Normal to mildly increased','Moderately increased','Severely increased'])}),hide_index=True,use_container_width=True)
        st.subheader(T('Matriz de categorías pronósticas','Prognostic category matrix'), anchor=False)
        # Original HTML rendering of category facts, not a copied KDIGO image.
        html='<div style="overflow-x:auto"><table style="width:100%;border-collapse:separate;border-spacing:4px;text-align:center;font-size:0.9rem"><thead><tr><th>G / A</th><th>A1</th><th>A2</th><th>A3</th></tr></thead><tbody>'
        for i,g in enumerate(G_NAMES):
            html+=f'<tr><th style="padding:8px">{g} · {G_RANGES[i].replace("<","&lt;")}</th>'
            for j,risk in enumerate(RISK[i]):
                chosen=(gi==i and ai==j); border='3px solid #193c55' if chosen else '3px solid transparent'
                mark=' · '+T('Resultado','Result') if chosen else ''
                html+=f'<td style="background:{COLORS[risk]};color:#20282e;padding:12px 6px;border:{border};border-radius:6px;font-weight:{700 if chosen else 500}">{risks[risk]}{mark}</td>'
            html+='</tr>'
        st.markdown(html+'</tbody></table></div>',unsafe_allow_html=True)
        if gi is not None and ai is not None: st.write(f'**{G_NAMES[gi]} / A{ai+1} — '+risks[RISK[gi][ai]]+'** · '+T('categorías de los valores ingresados, no diagnóstico ni probabilidad individual.','categories of entered values, not a diagnosis or individual probability.'))
        else: st.caption(T('Introduce eTFG y RAC válidas para resaltar una celda.','Enter valid eGFR inputs and ACR to highlight a cell.'))
        st.warning(T('La matriz se interpreta en el contexto de ERC. G1/G2 con A1, sin otros marcadores de daño renal, no constituyen ERC. No asumir cronicidad con una sola medición.','The matrix is interpreted in the context of CKD. G1/G2 with A1 and no other markers of kidney damage do not constitute CKD. Do not assume chronicity from a single measurement.'))
        st.caption(T('Fuente: KDIGO 2024, S126–S127 y tablas 1–3, S137. Límites de RAC por unidad según la guía (equivalencias aproximadas).','Source: KDIGO 2024, S126–S127 and Tables 1–3, S137. Unit-specific ACR thresholds follow the guideline (approximate equivalence).'))
    with tabs[4]:
        st.header(T('Calidad analítica y valores de referencia','Analytical quality and reference intervals'), anchor=False)
        st.dataframe(pd.DataFrame({T('Marcador','Marker'):T(['Creatinina','Cistatina C'],['Creatinine','Cystatin C']),T('Trazabilidad','Traceability'):['IDMS','ERM-DA471/IFCC'],'CV':['<2.3%','<2.0%'],T('Sesgo deseable','Desirable bias'):['<3.7%','<3.2%']}),hide_index=True,use_container_width=True)
        st.write(T('Objetivos analíticos de KDIGO, tabla 11, S184. Se prefiere creatinina enzimática cuando sea posible; medir ambos marcadores en la misma muestra. La aplicación no verifica calibración ni calidad de las muestras.','Analytical targets from KDIGO, Table 11, S184. Enzymatic creatinine is preferred where possible; measure both markers on the same sample. The app does not verify calibration or sample quality.'))
        st.subheader(T('Unidades e intervalos','Units and reference intervals'), anchor=False)
        st.dataframe(pd.DataFrame({T('Variable','Variable'):T(['Creatinina','Cistatina C','RAC','eTFG'],['Creatinine','Cystatin C','ACR','eGFR']),T('Unidades','Units'):['mg/dL / µmol/L','mg/L','mg/g / mg/mmol',UNITS],T('Referencia','Reference'):T(['Según laboratorio y método','Según laboratorio y método','Categorías A1–A3','Categorías G1–G5'],['Laboratory and method specific','Laboratory and method specific','A1–A3 categories','G1–G5 categories'])}),hide_index=True,use_container_width=True)
        st.write(T('Creatinina: µmol/L ÷ 88.4 = mg/dL. RAC: mg/g × 0.113 ≈ mg/mmol. No se asignan intervalos universales de creatinina o cistatina C sin documentación del método.','Creatinine: µmol/L ÷ 88.4 = mg/dL. ACR: mg/g × 0.113 ≈ mg/mmol. No universal creatinine or cystatin C intervals are assigned without method documentation.'))
        st.subheader(T('Limitaciones','Limitations'), anchor=False)
        st.write(T('La eTFG no es una TFG medida. Esta herramienta está dirigida a adultos con función renal relativamente estable; no se diseñó para embarazo, lesión renal aguda, diálisis o población pediátrica. En trasplante renal o enfermedad grave se requiere interpretación especializada. Masa muscular, dieta, amputaciones y ciertos fármacos afectan la creatinina; inflamación, obesidad, tabaquismo, alteraciones tiroideas y corticoides pueden afectar la cistatina C. No usar para ajustar fármacos de manera autónoma.','eGFR is not measured GFR. This tool is intended for adults with relatively stable kidney function; it is not designed for pregnancy, acute kidney injury, dialysis or children. Transplant recipients or severe illness require specialist interpretation. Muscle mass, diet, amputations and some drugs affect creatinine; inflammation, obesity, smoking, thyroid disorders and corticosteroids may affect cystatin C. Do not use for independent drug dosing.'))
    with tabs[5]:
        st.header(T('Fundamento científico','Scientific basis'), anchor=False)
        st.write(T('Inker 2012 desarrolló las ecuaciones con 5,352 participantes y validó externamente en 1,119. En esa validación, P30 fue 85.9% para cistatina C 2012 y 91.5% para la combinada 2012 (no la ecuación 2021). P30 indica la proporción de estimaciones dentro del 30% de la TFG medida, no la probabilidad de que un resultado individual sea correcto.','Inker 2012 developed equations with 5,352 participants and externally validated them in 1,119. In that validation, P30 was 85.9% for 2012 cystatin C and 91.5% for the 2012 combined equation (not the 2021 equation). P30 is the proportion within 30% of measured GFR, not the probability that an individual result is correct.'))
        st.write(T('Las ecuaciones CKD-EPI 2021 fueron desarrolladas sin coeficiente racial y la estimación combinada mostró mejor exactitud global que las estimaciones individuales en la validación original. KDIGO recomienda la ecuación combinada cuando está disponible cistatina C. El software no cuenta con validación clínica propia ni específica para México.','The 2021 CKD-EPI equations were developed without a race coefficient; the combined estimate had better overall accuracy than single-marker estimates in original validation. KDIGO recommends the combined equation when cystatin C is available. This software has no independent clinical or Mexico-specific validation.'))
        st.subheader(T('Referencias bibliográficas','References'), anchor=False)
        st.markdown('''1. Inker LA, et al. Estimating Glomerular Filtration Rate from Serum Creatinine and Cystatin C. *N Engl J Med*. 2012;367:20–29. [DOI](https://doi.org/10.1056/NEJMoa1114248).\n2. Inker LA, et al. New Creatinine- and Cystatin C–Based Equations to Estimate GFR without Race. *N Engl J Med*. 2021;385:1737–1749. [DOI](https://doi.org/10.1056/NEJMoa2102953).\n3. Miller WG, et al. National Kidney Foundation Laboratory Engagement Working Group Recommendations for Implementing the CKD-EPI 2021 Race-Free Equations for Estimated Glomerular Filtration Rate: Practical Guidance for Clinical Laboratories. *Clin Chem*. 2022;68:511–520. [DOI](https://doi.org/10.1093/clinchem/hvab278).\n4. KDIGO CKD Work Group. KDIGO 2024 Clinical Practice Guideline for the Evaluation and Management of Chronic Kidney Disease. *Kidney Int*. 2024;105(4S):S117–S314. [DOI](https://doi.org/10.1016/j.kint.2023.10.018).\n5. NIDDK. [eGFR Equations for Adults]('''+NIDDK+''').''')
        st.subheader(T('Privacidad y alcance','Privacy and scope'), anchor=False)
        st.write(T('La aplicación no implementa bases de datos, exportación ni registro de valores clínicos. Los datos se procesan en la sesión de Streamlit; el proveedor de alojamiento puede tratar información técnica. Evita datos identificables.','The app does not implement databases, exports or clinical-value logging. Inputs are processed in the Streamlit session; the hosting provider may process technical information. Avoid identifying data.'))
    st.divider()
    st.markdown('**'+T('Autores','Authors')+'**')
    st.write(' · '.join(AUTHORS))
    st.caption(T('Proyecto académico de medicina y química clínica. Versión preliminar ','Academic project in medicine and clinical chemistry. Preliminary version ')+VERSION+' · '+T('Última actualización: ','Last updated: ')+UPDATED)
    st.caption(T('No sustituye la valoración médica. No es un dispositivo médico validado.','Does not replace medical assessment. Not a validated medical device.'))


if __name__ == '__main__':
    main()
