using namespace std;

double DoubleUniformDistribution(double min, double max);

struct Transition
{
	int id;
	double time;
	Transition()
	{
		id = -1;
		time = 0.0;
	};
	Transition(int new_id, double newTime)
	{
		id = new_id;
		time = newTime;
	};
	void Update(double newTime) { time = newTime; };
};

class EventQueue
{
private:
	list<Transition> c;

public:
	EventQueue(list<Transition> newc) { c = newc; };
	EventQueue() { c = list<Transition>(); };
	Transition top() { return *c.begin(); };
	void pop() { c.pop_front(); };
	void push_back(Transition t) { c.push_back(t); };
	void erase(int id)
	{
		int n = c.size();
		for (list<Transition>::iterator it = c.begin(); it != c.end(); ++it)
			if (it->id == id)
			{
				c.erase(it);
				break;
			}
	};
	void erase_queue() { c = list<Transition>(); };
	bool empty() { return c.empty(); };
	int getSize() { return c.size(); };
	Transition last()
	{
		if (c.empty())
			return Transition();
		int n = c.size();
		list<Transition>::iterator it = c.begin();
		advance(it, n - 1);
		return *it;
	};
	~EventQueue() {}
	void push(Transition t)
	{
		if (c.empty())
		{
			c.push_front(t);
			return;
		}
		list<Transition>::iterator low = c.begin();
		list<Transition>::iterator high = c.end();
		list<Transition>::iterator mid;

		int d = distance(low, high);

		while (d > 1)
		{
			mid = next(low, d / 2);
			if (t.time < mid->time)
				high = mid;
			else
				low = mid;
			d = distance(low, high);
		}

		if (t.time < low->time)
			c.insert(low, t);
		else
			c.insert(high, t);
	};
};


struct Propensities
{
	int R; // number of reactions
	int n; // n is the dimension of the lattice
	vector<vector<double>> Reactions;
	Propensities()
	{
		R = 0;
		n = 0;
		Reactions = vector<vector<double>>();
	};
	Propensities(int newR, int newN)
	{
		R = newR;
		n = newN;
		Reactions = vector<vector<double>>(R, vector<double>(newN, 0.0));
	};
	void swap(Propensities &obj)
	{
		R = obj.R;
		n = obj.n;
		Reactions = obj.Reactions;
	};
	Propensities &operator=(Propensities a)
	{
		this->swap(a);
		return *this;
	};
	int wchEvent(int voxel)
	{
		double rand = DoubleUniformDistribution(0, 1), sum = 0.0, sumTot = 0.0;
		for (int i = 0; i < R; i++)
			sumTot += Reactions[i][voxel - 1];
		int reactionID = 0;
		while (reactionID < R)
		{
			sum += Reactions[reactionID++][voxel - 1];
			if (rand < sum / sumTot)
				return reactionID;
		}
		return -1; //error code
	};

	void clear(int voxel) // voxel indexing starts from 1
	{
		for (int i = 0; i < R; i++)
			Reactions[i][voxel - 1] = 0.0;
	};
};
