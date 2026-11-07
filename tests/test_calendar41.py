import unittest


class StaggeredCalendars(unittest.TestCase):
    def test_complete_instance_has_twelve_physical_fleets(self):
        from calendar41.mechanism import enumerate_instance
        r=enumerate_instance()
        self.assertEqual(r['candidates'],100)
        self.assertEqual(len(r['feasible']),12)
        self.assertEqual(len(r['normalized']),4)
        self.assertEqual(set(str(v['check']['base']['objective'])for v in r['feasible']),{'56'})
        self.assertTrue(all(v['check']['transport']==40 and v['check']['weighted_release']==62 for v in r['feasible']))

    def test_group_gaps_suffice_without_a_common_gap(self):
        from calendar41.mechanism import instance,calendars,free_gap
        self.assertEqual(free_gap(calendars()[0],6,15,8),5)
        self.assertEqual(free_gap(calendars()[1],6,15,8),5)
        self.assertEqual(free_gap(calendars()[0]+calendars()[1],6,15,8),3)

    def test_suffix_exchange_preserves_base_but_violates_both_calendars(self):
        from calendar41.mechanism import instance,calendars,fleet,validate
        from review88_physical.physical import validate_fleet
        p=fleet(10,8,old=(0,1),new=(3,2))
        self.assertTrue(validate_fleet(instance(),p)['valid'])
        r=validate(instance(),calendars(),p)
        self.assertFalse(r['valid'])
        self.assertEqual(set(r['calendar_conflicts']),{0,1})

    def test_missing_future_verification_has_no_verified_release(self):
        from calendar41.mechanism import instance,calendars,fleet,validate
        c=[v[:-1]for v in calendars()]
        r=validate(instance(),c,fleet(8,10))
        self.assertFalse(r['valid'])
        self.assertIn('missing subsequent verification',r['errors'])

    def test_every_feasible_fleet_normalizes_without_cost_or_service_changes(self):
        from calendar41.mechanism import enumerate_instance
        r=enumerate_instance()
        for v in r['feasible']:
            self.assertTrue(v['normalized_check']['valid'])
            self.assertEqual(v['check']['base']['objective'],v['normalized_check']['base']['objective'])
            self.assertEqual(v['check']['releases'],v['normalized_check']['releases'])
            for p,q in zip(v['plans'],v['normalized_plans']):
                self.assertEqual(p['routes'],q['routes']);self.assertEqual(p['assays'],q['assays'])


if __name__=='__main__':unittest.main()
