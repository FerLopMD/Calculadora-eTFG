import unittest
import math
from app import estimate, reported, g_index, a_index, parse_value, RISK

class TestCalculator(unittest.TestCase):
    def test_published_creatinine_table(self):
        # Miller et al. 2022, Table 2, p. 515; results rounded to integers.
        expected={'M':[[142,110,68,48],[118,92,56,40],[98,76,47,33]],'F':[[132,83,51,36],[109,69,42,30],[91,57,35,25]]}
        for sex,rows in expected.items():
            for age,row in zip((20,50,80),rows):
                for cr,want in zip((.6,1,1.5,2),row):
                    with self.subTest(sex=sex,age=age,cr=cr):
                        self.assertEqual(reported(estimate(age,sex,cr)['cr']['value']),want)
    def test_combined_piecewise(self):
        # Miller 2022 Table 4 independent piecewise expression.
        for sex,k,low,sf in [('F',.7,-.219,.963),('M',.9,-.144,1)]:
            for cr in (k-.01,k,k+.01,1.5):
                for cys in (.79,.8,.81,1.8):
                    with self.subTest(sex=sex,cr=cr,cys=cys):
                        want=135*(cr/k)**(low if cr<=k else -.544)*(cys/.8)**(-.323 if cys<=.8 else -.778)*.9961**55*sf
                        self.assertAlmostEqual(estimate(55,sex,cr,cys)['both']['value'],want,places=10)
    def test_cystatin_piecewise(self):
        # Inker 2012, Table 2.
        for sex,sf in [('F',.932),('M',1)]:
            for cys in (.79,.8,.81,1.8):
                want=133*(cys/.8)**(-.499 if cys<=.8 else -1.328)*.996**55*sf
                self.assertAlmostEqual(estimate(55,sex,cys=cys)['cys']['value'],want,places=10)
    def test_creatinine_knots(self):
        for sex,k,alpha,sf in [('M',.9,-.302,1),('F',.7,-.241,1.012)]:
            for cr in (k-.01,k,k+.01):
                want=142*(cr/k)**(alpha if cr<=k else -1.2)*.9938**55*sf
                self.assertAlmostEqual(estimate(55,sex,cr)['cr']['value'],want,places=10)
    def test_sample_and_units(self):
        a=estimate(55,'M',1.5,1.8)
        self.assertAlmostEqual(a['both']['value'],43.884560675757406)
        b=estimate(55,'M',132.6/88.4,1.8)
        self.assertAlmostEqual(a['both']['value'],b['both']['value'])
        self.assertEqual(set(estimate(55,'M',cr=1.5)),{'cr'})
        self.assertEqual(set(estimate(55,'M',cys=1.8)),{'cys'})
    def test_categories(self):
        for v,i in [(90,0),(89.99,1),(60,1),(59.99,2),(45,2),(44.99,3),(30,3),(29.99,4),(15,4),(14.99,5)]:
            self.assertEqual(g_index(v),i)
        for unit,lo,hi in [('mg/g',30,300),('mg/mmol',3,30)]:
            for v,i in [(0,0),(lo-.01,0),(lo,1),(hi,1),(hi+.01,2)]: self.assertEqual(a_index(v,unit),i)
        self.assertIsNone(a_index(None,'mg/g'))
        self.assertEqual(RISK,((0,1,2),(0,1,2),(1,2,3),(2,3,3),(3,3,3),(3,3,3)))
        self.assertEqual(reported(59.5),60)
        self.assertEqual(reported(59.49),59)
    def test_invalid(self):
        for x in ('nan','inf','-1','0','abc',''):
            with self.assertRaises(ValueError): parse_value(x,'x')
        self.assertEqual(parse_value('1,5','x'),1.5)
        self.assertIsNone(parse_value('','x',True))
        self.assertEqual(parse_value('0','x',True,True),0)
        for args in [(17,'M',1,1),(55,'',1,1),(55,'F',None,None),(55,'M',0,1),(55,'F',1,float('nan'))]:
            with self.assertRaises(ValueError): estimate(*args)

if __name__=='__main__': unittest.main()
