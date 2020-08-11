using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MPS_Method : MonoBehaviour
{
    //モデル粒子の情報を呼び出す
    public GameObject particle;

    //水の密度(kg/m^3)を保存するリスト(5~30℃)
    public static float[] densities = new float[26] { 999.993f, 999.974f, 999.938f, 999.987f, 999.887f, 999.821f, 999.741f, 999.647f, 999.539f, 999.418f, 999.284f, 999.138f, 998.980f, 998.628f, 998.436f, 998.223f, 998.019f, 997.794f, 997.560f, 997.316f, 997.062f, 996.799f, 996.526f, 996.244f, 995.954f, 995.654f };

    //水の粘度(Pa・s)を保存するリスト(5~30℃)
    public static float[] viscosities = new float[26] { 1518.94f, 1472.21f, 1427.74f, 1385.39f, 1345.02f, 1306.51f, 1269.74f, 1234.61f, 1201.02f, 1168.87f, 1138.09f, 1108.59f, 1080.30f, 1053.16f, 1027.11f, 1002.07f, 978.01f, 954.86f, 932.59f, 911.15f, 890.49f, 870.57f, 851.37f, 832.85f, 814.96f, 797.70f };

    //水温を標準入力から受け取る(5~30℃)
    public int temperature;

    //影響半径を標準入力から受け取る
    public float r_e;

    //位置を保存するリスト,速度を保存するリスト,粒子数密度を保存するリストの定義
    public static List<Vector3> position_l = new List<Vector3>();
    public static List<Vector3> velocity_l = new List<Vector3>();
    public static List<float> n_l = new List<float>();

    //流体を表す粒子数の数
    public static int cnt = 0;
    //床を表す粒子の数
    public static int add_cnt = 0;

    //λ
    public static float lambda = 0f;

    //粒子数密度の初期値
    public static float n0 = 0f;

    //重み関数
    float W(float r)
    {
        if (r_e < r) return 0f;
        return r_e / r - 1;
    }
    //二乗を計算する関数
    float Pow2(float x)
    {
        return x * x;
    }

    //最初一回だけ呼び出される関数
    void Start()
    {
        //等間隔で粒子を配置していく(隙間はなし、体心立方格子みたいな感じ)
        //(x,y,z)で粒子の中心座標を示す
        //半径0.05f
        //xの範囲
        for (var x = -2; x < 3; x++)
        {
            //yの範囲
            for (var y = 0; y < 4; y++)
            {
                //zの範囲
                for (var z = -2; z < 3; z++)
                {
                    //モデル粒子を(x,y,z)の点に回転なしで作成
                    Instantiate(particle, new Vector3(0.1f * x, 0.1f * y, 0.1f * z), Quaternion.identity);

                    //それぞれのリストの初期化
                    Vector3 pos = new Vector3(0.1f * x, 0.1f * y, 0.1f * z);
                    Vector3 vel = new Vector3(0f, 0f, 0f);
                    float num = 0f;
                    position_l.Add(pos);
                    velocity_l.Add(vel);
                    n_l.Add(num);

                    //粒子数のカウント
                    cnt++;
                }
            }
        }

        //粒子で床を作ろう
        //等間隔で平たく作っていく
        //こいつらは動かさないので計算量の心配はせずに粒子数を増やせる
        //厚さは影響半径よりも大きくなるようにする
        //上と同様の考え方で
        for (var x = -30; x < 31; x++)
        {
            for (var y = -2; y < 0; y++)
            {
                for (var z = -30; z < 31; z++)
                {
                    Instantiate(particle, new Vector3(0.1f * x, 0.1f * y, 0.1f * z), Quaternion.identity);

                    //床粒子のもきちんと加えておこう
                    Vector3 pos = new Vector3(0.1f * x, 0.1f * y, 0.1f * z);
                    Vector3 vel = new Vector3(0f, 0f, 0f);
                    float num = 0f;
                    position_l.Add(pos);
                    velocity_l.Add(vel);
                    n_l.Add(num);

                    //粒子数のカウント
                    add_cnt++;
                }
            }
        }

        //λは使いまわせるからここで計算してしまおう
        //λを表す分数の分子と分母の定義(λ=n/d)
        float n = 0f, d = 0f;

        //Σi≠j(i=52)
        float xi_x = position_l[52].x;
        float xi_y = position_l[52].y;
        float xi_z = position_l[52].z;
        for (int j = 0; j < cnt + add_cnt; j++)
        {
            //粒子xjの座標の取得
            if (52 == j) continue;
            float xj_x = position_l[j].x;
            float xj_y = position_l[j].y;
            float xj_z = position_l[j].z;

            //n0に関する計算
            n0 += W((float)Math.Sqrt(Pow2(xj_x - xi_x) + Pow2(xj_y - xi_y) + Pow2(xj_z - xi_z)));

            //λに関する計算
            //分子nについての計算
            n += (Pow2(xj_x - xi_x) + Pow2(xj_y - xi_y) + Pow2(xj_z - xi_z)) * W((float)Math.Sqrt(Pow2(xj_x - xi_x) + Pow2(xj_y - xi_y) + Pow2(xj_z - xi_z)));

            //分母dについての計算
            d += W((float)Math.Sqrt(Pow2(xj_x - xi_x) + Pow2(xj_y - xi_y) + Pow2(xj_z - xi_z)));
        }
        //λの算出
        lambda = n / d;
    }

    //物理演算
    //0.02秒毎に呼び出す
    void FixedUpdate()
    {
        //密度と粘度の取得
        float density = densities[temperature - 5];
        float viscosity = viscosities[temperature - 5];

        //各粒子に対して計算
        for (int i = 0; i < cnt; i++)
        {
            //粘性項を離散化したときに出てくるΣ
            //各軸ごとに計算するので3つ
            float sigma_x = 0f;
            float sigma_y = 0f;
            float sigma_z = 0f;

            //粒子xiの座標の取得
            float xi_x = position_l[i].x;
            float xi_y = position_l[i].y;
            float xi_z = position_l[i].z;

            //粒子xiの速度の取得
            float xi_vx = velocity_l[i].x;
            float xi_vy = velocity_l[i].y;
            float xi_vz = velocity_l[i].z;

            //Σi≠j
            for (int j = 0; j < cnt + add_cnt; j++)
            {
                //粒子xjの座標の取得
                if (i == j) continue;
                float xj_x = position_l[j].x;
                float xj_y = position_l[j].y;
                float xj_z = position_l[j].z;

                //粒子xjの速度の取得
                float xj_vx = velocity_l[j].x;
                float xj_vy = velocity_l[j].y;
                float xj_vz = velocity_l[j].z;

                //Σの計算
                sigma_x += (xj_vx - xi_vx) * W((float)Math.Sqrt(Pow2(xj_x - xi_x) + Pow2(xj_y - xi_y) + Pow2(xj_z - xi_z)));
                sigma_y = (xj_vy - xi_vy) * W((float)Math.Sqrt(Pow2(xj_x - xi_x) + Pow2(xj_y - xi_y) + Pow2(xj_z - xi_z)));
                sigma_z = (xj_vz - xi_vz) * W((float)Math.Sqrt(Pow2(xj_x - xi_x) + Pow2(xj_y - xi_y) + Pow2(xj_z - xi_z)));
            }

            //粘性項の算出(d = 3)
            float viscosity_x = 2 * 3 * sigma_x / lambda / n0;
            float viscosity_y = 2 * 3 * sigma_y / lambda / n0;
            float viscosity_z = 2 * 3 * sigma_z / lambda / n0;

            //速度と位置の更新(仮)
            //漸化式に倣って更新していく(Δt = 0.02, g = -9.81)
            //速度の更新(仮)
            float vx_temporary = velocity_l[i].x + 0.02f * (viscosity * viscosity_x / density - 9.81f);
            float vy_temporary = velocity_l[i].y + 0.02f * (viscosity * viscosity_y / density - 9.81f);
            float vz_temporary = velocity_l[i].z + 0.02f * (viscosity * viscosity_z / density - 9.81f);
            velocity_l[i] = new Vector3(vx_temporary, vy_temporary, viscosity_z);

            //位置の更新(仮)
            float x_temporary = position_l[i].x + 0.02f * vx_temporary;
            float y_temporary = position_l[i].y + 0.02f * vy_temporary;
            float z_temporary = position_l[i].z + 0.02f * vz_temporary;
            position_l[i] = new Vector3(x_temporary, y_temporary, z_temporary);
        }

        //仮位置における粒子数密度nの算出
        //仮位置における粒子数密度n
        float n = 0f;

        //Σi≠j
        for (int i = 0; i < cnt; i++)
        {
            //粒子xiの座標の取得
            float xi_x = position_l[i].x;
            float xi_y = position_l[i].y;
            float xi_z = position_l[i].z;

            for (int j = 0; j < cnt + add_cnt; j++)
            {
                //粒子xjの座標の取得
                if (i == j) continue;
                float xj_x = position_l[j].x;
                float xj_y = position_l[j].y;
                float xj_z = position_l[j].z;

                //nに関する計算
                n += W((float)Math.Sqrt(Pow2(xj_x - xi_x) + Pow2(xj_y - xi_y) + Pow2(xj_z - xi_z)));
            }
        }
    }
}
