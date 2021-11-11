#pragma once

class B4P_In {
public:
	virtual void Nothing(void)=0;	// 抽象クラスであることを示すため，意味のない関数を定義する．
	
public:
	int ballnum;					// 転動体の個数
	double balldia;				// 基本玉径（初期値計算のみに用いる）
	double ballpcd;				// 転動体pcd（『4つの溝R中心pcdの平均値』を取る．初期値計算のみに使用）
	double cos_alp0;			// 接触角の余弦（初期値計算のみに使用）

	struct Tribology {
		enum RollingResistance {
			RollingResistanceNothing = 0, // 転がり粘性抵抗なし
			GoksemAihara = 1,	// Goksem-相原の式（baltac実装式）
			Aihara = 2,			// 相原の式
			Fujiwara = 3,		// 藤原の式
			Houpert = 4,		// Houpertの式
		} rollingresistance;

		enum Coulomb {
			CoulombNothing = 0,	// クーロン摩擦なし．
			Tangent = 1,	// タンジェントカーブ．
		} coulomb;
		double coulomb_slope;

		enum FilmThickness {
			FilmThicknessNothing = 0,	// 油膜なし．
			HamrockDowsonHc = 1,			// HamrockDowson式の中央油膜厚さを摩擦計算に使用．
			HamrockDowsonHmin = 2,			// HamrockDowson式の最小油膜厚さを摩擦計算に使用．
		} filmThickness;

		enum Hysteresis {
			HysteresisNothing = 0,		// ヒステリシスなし．
			Kakuta = 1,					// 角田の式
		}hysteresis;
		double hysteresis_factor;		// ヒステリシス損失係数
	} TB;

	enum CageType {
		spherical_machined_cage		= 0,	// 球面もみ抜き
		snap_cage					= 1,	// 冠型
		cylindrical_machined_cage	= 2,	// 円筒もみ抜き
		corrugated_press_cage		= 3,	// 波型プレス
		cylindrical_snap_cage		= 4,	// 円筒冠型TY
		snap_press_cage				= 5,	// 冠型プレス
		full_ball					= -2,	// 総玉軸受
		no_cage						= -1	// 保持器無し
	} cage_type;
	   
	// 転動体の密度，剛性，玉径などは玉番号ごとに定義
	struct Ball {
		double den;				// 転動体密度[kg/m^3]
		double E;				// 転動体剛性[Pa]
		double por;				// 転動体ポアソン比
		double rms;				// 転動体粗さrms
		double dia;				// 転動体外径（転動体毎に定義）
	} *BL;
									   

									   
	struct Ring {
		double den;				// 内外輪密度[kg/m^3]
		double E;				// 内外輪ヤング率[Pa]
		double por;				// 内外輪ポアソン比
		double rms;				// 内外輪粗さrms
		double R[2];			// 内外輪溝R半径[m] 0:-x側溝, 1:+x側溝
		double Rox[2];			// 内外輪溝R中心x座標[m]
		double Rod[2];			// 内外輪溝R中心直径[m]
		double hedge[2];		// 内外輪溝肩高さ[m]
		double m;				// 内外輪重量[kg]
		double Ix;				// 内外輪慣性モーメントx方向[kg・m^2]
		double Iyz;				// 内外輪慣性モーメントy方向[kg・m^2]
	} IR, OR;
	
	struct BallRingPair {
		double mu;				// 摩擦係数	（転動体-内輪間）
		double dzeta;			// 減衰比		（転動体-内輪間）
	} BIP, BOP;

	struct BallCagePair {
		double mu;				// 摩擦係数	（転動体-保持器間）
		double dzeta;			// 減衰比		（転動体-保持器間）
	} BCP;

	struct CageRingPair {
		double mu;				// 摩擦係数	（保持器-内輪間）
		double dzeta;			// 減衰比		（保持器-内輪間）
	}CIP, COP;

	// 冠型保持器
	struct SnapCage {
		double dout;			// 冠型保持器外径[m]			// initで使う
		double din;				// 冠型保持器内径[m]			// initで使う
		double h;				// 冠型保持器高さ[m]			// initで使う
		double jc;				// 冠型保持器jc量[m] ※（保持器のPCD）－（軸受のPCD）	// initで使う
		double R;				// 冠型保持器ポケット直径隙間[m](pocket clearance)					// initで使う
		double ropen;			// 冠型保持器ポケット開口部半径[m](pocket opening Radius)			// initで使う
		double Kface;			// 転動体-保持器剛性（面接触）
		double Kedgein;			// 転動体-保持器剛性（エッジ内側）
		double Kedgeout;		// 転動体-保持器剛性（エッジ外側）
		double Kcornerin;		// 転動体-保持器剛性（角内側）
		double Kcornerout;		// 転動体-保持器剛性（角外側）
		double Kopen;			// 転動体-保持器剛性（開口部）
	} Snap;

	struct Cage {
		double den;				// 冠型保持器密度[kg/m^3]	// 接触剛性自動計算以外で使わない
		double E;				// 冠型保持器ヤング率[Pa]	// 接触剛性自動計算以外で使わない
		double por;				// 冠型保持器ポアソン比		// 接触剛性自動計算以外で使わない
		double rms;				// 冠型保持器粗さrms			// 摩擦係数自動計算以外で使わない
		double rmg0[3];			// 冠型保持器　重心を基準にした幾何中心の位置ベクトル（Y・Z成分を0にしたもの）
		double m;				// 冠型保持器重量[kg]
		double Ix;				// 冠型保持器慣性モーメントx方向[kg・m^2]
		double Iyz;				// 冠型保持器慣性モーメントy方向・z方向[kg・m^2]
	} Cage;
			
	struct Lubrication {
		double eta0;			// 粘度 [Pa*s]
		double beta0;			// 温度粘度係数 [K^-1]
		double k0;				// 油熱伝導率 [W/(K*m)]
		double alpha0;			// 圧力粘度係数 [Pa^-1]
		double lm0;				// メニスカス長さ [m]
	} LB;

	double LoadIn[6];			// 内輪荷重x,y,z,Rx,Ry,Rz[N][N・m]
		   
	double omegair;				// 内輪回転速度[rad/s]
	double omegaor;				// 外輪回転速度[rad/s]
	int msmax = 21;				// 接触楕円分割数
	
	struct Rigid {
		double l;
		double t;
		double g[3];			// 重力(x,y,z)
	} rigid;

	struct Bound {
		bool v_const[3];
		bool w_const[3];
	} bound;
};





































	// bool autocalc_Inner_m;		// 内輪重量自動計算（true:自動計算，false:自動計算しない）
	// bool autocalc_Inner_Ix;		// 内輪慣性モーメントx方向自動計算（true:自動計算，false:自動計算しない）
	// bool autocalc_Inner_Iy;		// 内輪慣性モーメントy方向自動計算（true:自動計算，false:自動計算しない）
	// bool gravity_pulls_ring;	// 軌道輪への重力考慮(true:重力考慮する，false:重力考慮しない，デフォルト：true)
	// bool autocalc_Snap_m;		// 冠型保持器　重量自動計算（true:自動計算，false:自動計算しない）
	// bool autocalc_Snap_Ix;		// 冠型保持器　慣性モーメントx方向自動計算（true:自動計算，false:自動計算しない）
	// bool autocalc_Snap_IyIz;	// 冠型保持器　慣性モーメントy方向・z方向自動計算（true:自動計算，false:自動計算しない）
	// Vector3d rmg_;			// 保持器の重心から見た幾何中心のy方向位置[m]?
	// Vector3d rmg0_;			// 保持器の重心から見た幾何中心のy方向位置[m]
	// bool autocalc_Snap_rmg[3];	// 冠型保持器　幾何中心自動計算（true:自動計算，false:自動計算しない）
	// DParam Outer;
	// bool autocalc_Outer_m;		// 外輪重量自動計算（true:自動計算，false:自動計算しない）
	// bool autocalc_Outer_Ix;		// 外輪慣性モーメントx方向自動計算（true:自動計算，false:自動計算しない）
	// bool autocalc_Outer_Iy;		// 外輪慣性モーメントy方向自動計算（true:自動計算，false:自動計算しない）
	// DParam ContaBI;
	// DParam ContaBO;
	// DParam ContaBC;
	// double BC_Kinput;		// 剛性入力値	（転動体-保持器間）
	// DParam ContaCI;
	// DParam ContaCO;	





