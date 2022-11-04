def change_age(data):
    idade = [i for i in data.columns if('AGE' in i)]
    idades = [list(i) for i in data[idade].fillna(-1).values]
    idades = [i[:i.index(-1)] if(-1 in i) else i for i in idades]
    for i in idade:
        data.drop(i, inplace=True, axis=1)
    data['AGE PERSONNE'] = idades
    return data
def change_mode(data):
    modes = [i for i in data.columns if('mode(s)' in i)]
    for i in range(len(modes)):
        mode = list(data.columns).index(modes[i])
        mode = data.columns[mode:mode+3]
        modes_ = [list(i) for i in data[mode].fillna(-1).values]
        modes_ = [i[:i.index(-1)] if(-1 in i) else i for i in modes_]
        for j in mode:
            data.drop(j, inplace=True, axis=1)
        data[modes[i]] = modes_
    return data 
def change_traches(data):
    tranches = [i for i in data.columns if('tranche(s)' in i)][0]
    tranche = list(data.columns).index(tranches)
    tranche_ = data.columns[tranche:tranche+5]
    tranche = [list(i) for i in data[tranche_].fillna(0).values]
    for i in tranche_:
        data.drop(i, inplace=True, axis=1)
    data[tranches] = tranche
    return data
def change_column(data):
    a = '''onda
    pergunta de número
    CONTROL Num Questionio versus recrutar Tel
    Tipo de questionário
    Idade do sujeito da pesquisa
    Sexo do sujeito da pesquisa
    Relação do entrevistado com a criança <15 anos Assunto da pesquisa
    Idade do entrevistado para a criança <15 anos Assunto da pesquisa
    Sexo do entrevistado para a criança <15 anos Sujeito da pesquisa
    Número de pessoas em casa
    Departamento
    Mais elevado grau
    situação profissional
    Em qual indústria você trabalha
    A criança vai à escola?
    criança cuidada em casa/com a família?
    criança mantida em assistência. Mat.\ n babá?
    nb filhos na assistente
    o assistente acolhe crianças em idade escolar?
    mantidos em um berçário?
    quantas crianças na creche
    frequência de creche
    Número de crianças em sua classe
    criança que come na cantina?
    ir para o acampamento diurno?
    ir para o acampamento de verão durante a escola?
    indo para o acampamento de verão durante as férias?
    profissão que leva a muitos contatos?
    número de contatos profissionais
    Tem ou não tem mais de 20 contactos profissionais
    Número de alunos da turma
    Estudante que come na cantina?
    dia1mês1
    Validação de CONTROLE dia1mês1 no período de onda 1 ou 2
    diasemana1
    Férias dia1mês1
    Número de contatos inseridos DIA 1
    PESSOA DE IDADE
    modo(s) de viagem semana
    modo(s) de viagem nós
    faixa(s) etária(s) de contatos no ambiente profissional
    '''
    colunas = a.split('\n')
    colunas = [i.replace('    ','') for i in colunas]
    data.columns = colunas[:-1]
    return data