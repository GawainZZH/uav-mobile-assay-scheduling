"""Finite stock-flow master; all rows keep exact rational coefficients."""
from fractions import Fraction as F
from types import SimpleNamespace as NS
from collections import defaultdict


def row(entries,rhs=0,eq=False,kind='linear'):
    return dict(entries={int(j):F(str(c))for j,c in entries.items()if c},rhs=F(str(rhs)),eq=bool(eq),kind=kind)


def add(entries,j,c):
    entries[j]=entries.get(j,F(0))+F(str(c))
    if not entries[j]:del entries[j]


def band_rows(m,i):
    b=m.bands[i];Fcols=m.family_columns[i]
    if b['infeasible']:
        rows=[row({j:1 for j in Fcols},0,False,'family_infeasible')]
        if i in m.z:rows.append(row({m.z[i]:1},0,False,'infeasible_family_zero_cost'))
        return rows
    if i not in m.z:return[]
    low={j:b['lower']for j in Fcols};add(low,m.z[i],-1)
    up={j:-b['upper']for j in Fcols};add(up,m.z[i],1)
    return[row(low,kind='family_lower'),row(up,kind='family_upper')]


def aggregate(m,form):
    e={}
    for kind in ('counts','u','theta','cost','boundary'):
        for key,c in form.get(kind,{}).items():
            if kind=='counts':
                if key not in m.count_columns:raise ValueError('unknown aggregate type')
                for j in m.count_columns[key]:add(e,j,c)
            elif kind=='u':
                if key not in m.u:raise ValueError('unknown threshold')
                add(e,m.u[key],c)
            elif kind in('theta','cost'):
                if key not in m.theta:raise ValueError('unknown stage cost')
                add(e,m.theta[key],c)
            else:
                if key not in m.boundary_columns:raise ValueError('unknown boundary flow')
                for j in m.boundary_columns[key]:add(e,j,c)
    return row(e,form['rhs'],form['eq'],'shared_aggregate')


def build(d,*,cutoff=None,forms=()):
    m=NS(domain=d,cost=[],upper=[],integer=[],rows=[],flow_columns=[],family_columns=defaultdict(list),
        count_columns=defaultdict(list),boundary_columns=defaultdict(list),z={},theta={},u={},bands=[dict(f)for f in d.families])
    def variable(cost,upper,integer=False):
        j=len(m.cost);m.cost.append(F(str(cost)));m.upper.append(F(str(upper)));m.integer.append(bool(integer));return j
    balance=defaultdict(dict);covers={j:{}for s in d.tasks for j in s};K=len(d.tasks)
    for a in d.arcs:
        i=a['family'];f=d.families[i];j=variable(0,1 if f['subset']else d.groups,True);m.flow_columns.append(j)
        m.family_columns[i].append(j);m.count_columns[f['stage'],(f['start'],f['end'],len(f['subset']))].append(j)
        add(balance[a['tail']],j,1);add(balance[a['head']]if a['head'][0]<K else balance['sink'],j,-1)
        if a['tail'][0]==0:m.boundary_columns[0,a['tail'][1]].append(j)
        m.boundary_columns[a['head'][0],a['head'][1]].append(j)
        for task in f['subset']:covers[task][j]=1
    source=(0,0,d.budget)
    if source not in balance:balance[source]={}
    for node,e in balance.items():m.rows.append(row(e,-d.groups if node=='sink'else d.groups if node==source else 0,True,'stock_flow'))
    for e in covers.values():m.rows.append(row(e,1,True,'cover'))
    for i,f in enumerate(d.families):
        if f['subset']:m.z[i]=variable(0,f['upper'])
    for k in range(K):m.theta[k]=variable(1,d.groups*max((f['upper']for f in d.families if f['stage']==k),default=0))
    for key in sorted(m.count_columns):
        for r in range(1,d.groups+1):m.u[key[0],key[1],r]=variable(0,1)
        e={j:-1 for j in m.count_columns[key]}
        for r in range(1,d.groups+1):add(e,m.u[key[0],key[1],r],1)
        m.rows.append(row(e,0,True,'type_count'))
        for r in range(1,d.groups):m.rows.append(row({m.u[key[0],key[1],r]:-1,m.u[key[0],key[1],r+1]:1},kind='threshold_order'))
    for k in range(K):
        e={m.theta[k]:-1}
        for i,f in enumerate(d.families):
            if f['stage']!=k:continue
            if f['subset']:add(e,m.z[i],1)
            elif not f['infeasible']:
                if not f['feasible']or f['lower']!=f['upper']:raise ValueError('empty family requires exact road plan')
                for j in m.family_columns[i]:add(e,j,f['lower'])
        m.rows.append(row(e,0,True,'stage_cost'))
    for i in range(len(d.families)):m.rows.extend(band_rows(m,i))
    for f in forms:m.rows.append(aggregate(m,f))
    if cutoff is not None:m.rows.append(row({j:1 for j in m.theta.values()},cutoff,False,'total_cutoff'))
    return m


def sparse(m):
    from scipy.sparse import coo_matrix
    rs=[];cs=[];vs=[]
    for i,r in enumerate(m.rows):
        for j,v in r['entries'].items():rs.append(i);cs.append(j);vs.append(float(v))
    return coo_matrix((vs,(rs,cs)),shape=(len(m.rows),len(m.cost))).tocsc()


def activity(r,x):return sum(float(c)*x[j]for j,c in r['entries'].items())-float(r['rhs'])
