#ifndef runMetadata_hh
#define runMetadata_hh 1

#include "TObject.h"

class runMetadata : public TObject
{
public:
	runMetadata() {;};
	virtual ~runMetadata() {;};

	inline void SetTotalEvents(int n) {totalEvents = n;};
	inline int GetTotalEvents() {return totalEvents;};

private:
	static runMetadata *theRunMetadata;

	int totalEvents;

	ClassDef(runMetadata,1);
};

#endif
